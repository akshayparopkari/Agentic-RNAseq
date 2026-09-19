import re
import json
import zipfile
import resource
import subprocess
import urllib.parse
from typing import Optional
from pathlib import Path

import yaml
import pandas as pd

_STAGE_SPECS = {
    "trim": {
        "output": lambda results_dir, sample, attempt, **_: f"{results_dir}/trim/{sample}_attempt{attempt}_trimmed.fastq",
        "requires": ("sample", "attempt"),
        "needs_attempt_config": False,  # attempt is a wildcard in the path itself
    },
    "fastqc": {
        "output": lambda results_dir, sample, attempt, **_: f"{results_dir}/qc/{sample}_attempt{attempt}_trimmed_fastqc.zip",
        "requires": ("sample", "attempt"),
        "needs_attempt_config": False,
    },
    "promote_trimmed": {
        "output": lambda results_dir, sample, **_: f"{results_dir}/trim/canonical/{sample}_trimmed.fastq",
        "requires": ("sample", "attempt"),
        "needs_attempt_config": True,  # accepted_attempt_path reads {sample}_attempt from --config
    },
    "align_count": {
        "output": lambda results_dir, sample, **_: f"{results_dir}/star/{sample}/{sample}ReadsPerGene.out.tab",
        "requires": ("sample", "attempt"),
        "needs_attempt_config": True,  # transitively depends on promote_trimmed
    },
    "collect_counts": {
        "output": lambda results_dir, **_: f"{results_dir}/counts/gene_raw_counts.txt",
        "requires": ("included_samples",),
        "needs_attempt_config": False,
    },
    "filter_metadata": {
        "output": lambda results_dir, **_: f"{results_dir}/metadata_filtered.xlsx",
        "requires": (),
        "needs_attempt_config": False,
    },
}


def _resolve_run_accession(identifier: str, metadata_xlsx: str) -> str:
    """
    Every pipeline-stage rule keys off the raw FASTQ filename (find_raw_fastq),
    not the metadata Sample_ID, so this has to be resolved before any target path is built.
    """
    df = pd.read_excel(metadata_xlsx)
    stems = df["FASTQ_file"].str.replace(r"\.fastq(\.gz)?$", "", regex=True)
    hit = df.index[
        (df["Sample_ID"] == identifier)
        | (stems == identifier)
        | stems.str.startswith(identifier)
    ]
    if len(hit) != 1:
        raise ValueError(
            f"Could not uniquely resolve '{identifier}' in {metadata_xlsx}; found {len(hit)} matches."
        )
    return stems[hit[0]]


def run_pipeline_stage(
    stage: str,
    sample: Optional[str] = None,
    attempt: Optional[int] = None,
    included_samples: Optional[list[str]] = None,
    reference_level: Optional[str] = None,
    configfile: str = "config.yaml",
    dry_run: bool = False,
    cores: int = 4,
    timeout: Optional[int] = None,
) -> dict:
    """
    Runs one named pipeline stage by deriving its real Snakemake target
    path(s) and required --config values from config.yaml and the
    Snakefile's own output conventions, instead of asking the caller to
    construct a raw target string from memory. An unknown stage or a
    missing required argument is a real error, not a guess.
    """
    if stage not in _STAGE_SPECS:
        return {
            "error": f"Unknown stage '{stage}'. Valid stages: {sorted(_STAGE_SPECS)}"
        }
    spec = _STAGE_SPECS[stage]

    provided = {
        "sample": sample,
        "attempt": attempt,
        "included_samples": included_samples,
        "reference_level": reference_level,
    }
    missing = [k for k in spec["requires"] if not provided.get(k)]
    if missing:
        return {
            "error": f"Stage '{stage}' requires {spec['requires']}; missing {missing}"
        }

    with open(configfile) as fh:
        cfg = yaml.safe_load(fh)
    results_dir = cfg["results_dir"]

    if sample is not None and stage in (
        "trim",
        "fastqc",
        "promote_trimmed",
        "align_count",
    ):
        try:
            sample = _resolve_run_accession(sample, cfg["metadata"])
        except ValueError as e:
            return {"error": str(e)}
    if included_samples and stage == "collect_counts":
        try:
            included_samples = [
                _resolve_run_accession(s, cfg["metadata"]) for s in included_samples
            ]
        except ValueError as e:
            return {"error": str(e)}
    targets = spec["output"](results_dir, sample=sample, attempt=attempt)
    targets = [targets] if isinstance(targets, str) else targets

    run_config = {}
    if spec["needs_attempt_config"]:
        run_config[f"{sample}_attempt"] = attempt
    if stage == "collect_counts":
        run_config["included_samples"] = ",".join(included_samples)
    if stage == "deseq2":
        run_config["reference_level"] = reference_level

    return run_snakemake_target(
        target=targets,
        config=run_config or None,
        configfile=configfile,
        dry_run=dry_run,
        cores=cores,
        timeout=timeout,
    )


def get_metadata_conditions(metadata_xlsx: str) -> dict:
    """Returns the distinct values in the Condition column of the metadata file."""
    import pandas as pd  # local import: only needed here

    df = pd.read_excel(metadata_xlsx)
    if "Condition" not in df.columns:
        return {"error": "metadata file has no 'Condition' column"}
    if len(df) == 0:
        return {"error": "metadata file has zero data rows"}
    return {
        "condition_levels": sorted(df["Condition"].dropna().unique().tolist()),
        "sample_count": len(df),
        "sample_ids": df.get("Sample_ID", df.get("SampleID")).tolist(),
    }


def resolve_sample_id(fastq_or_srr: str, metadata_xlsx: str) -> dict:
    df = pd.read_excel(metadata_xlsx)
    hits = df[
        (df["Sample_ID"].astype(str) == fastq_or_srr)
        | df["FASTQ_file"].str.startswith(fastq_or_srr)
    ]
    if len(hits) != 1:
        return {
            "error": f"Expected exactly one match for '{fastq_or_srr}', found {len(hits)}."
        }
    row = hits.iloc[0]
    return {
        "sample_id": row["Sample_ID"],
        "condition": row["Condition"],
        "fastq_file": row["FASTQ_file"],
    }


def parse_fastqc_summary(qc_dir: str, trimmed_fastq: str) -> dict:
    """
    Extract the metrics decide_qc_action needs, from FastQC's own output.
    Returns raw PASS/WARN/FAIL per module plus a couple of numeric proxies
    from fastqc_data.txt. Deliberately does NOT decide retry/exclude/proceed
    -- that's decide_qc_action's job, not this function's.

    trimmed_fastq must be the attempt-numbered path (e.g.
    results/trim/SRR123_attempt1_trimmed.fastq), not the canonical
    promoted one -- the fastqc rule keys off the attempt-numbered file.
    """
    stem = Path(trimmed_fastq).stem
    fastqc_zip = Path(qc_dir) / f"{stem}_fastqc.zip"
    if not fastqc_zip.exists():
        return {
            "error": f"expected {fastqc_zip} not found; did the fastqc target run first?"
        }

    with zipfile.ZipFile(fastqc_zip) as zf:
        data_txt = zf.read(f"{stem}_fastqc/fastqc_data.txt").decode(
            "utf-8", errors="replace"
        )
        summary_txt = zf.read(f"{stem}_fastqc/summary.txt").decode(
            "utf-8", errors="replace"
        )

    modules = {}
    for line in summary_txt.strip().splitlines():
        status, module, _ = line.split("\t")
        modules[module] = status  # "PASS" | "WARN" | "FAIL"

    total_sequences = _extract_field(data_txt, r"Total Sequences\t(\d+)")
    if total_sequences is None:
        raise ValueError(
            f"'Total Sequences' not found in {fastqc_zip}'s fastqc_data.txt -- "
            "can't make a depth call without a real read count. Defaulting to "
        )
    pct_gc = _extract_field(data_txt, r"%GC\t(\d+)")
    sequence_length = _extract_field(data_txt, r"Sequence length\t([\d\-]+)")

    return {
        "modules": modules,
        "total_sequences": int(total_sequences),
        "pct_gc": int(pct_gc) if pct_gc else None,
        "sequence_length": sequence_length,
        "adapter_content_status": modules.get("Adapter Content", "UNKNOWN"),
        "duplication_status": modules.get("Sequence Duplication Levels", "UNKNOWN"),
        "per_base_quality_status": modules.get("Per base sequence quality", "UNKNOWN"),
    }


def _extract_field(text: str, pattern: str) -> Optional[str]:
    m = re.search(pattern, text)
    return m.group(1) if m else None


def _ensure_file_descriptor_limit(min_soft: int = 4096) -> None:
    """
    STAR's BAM sorting opens threads * outBAMsortingBinsN temp files at
    once. Sets ulimit in process to account for this issue.
    """
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    target = min(min_soft, hard) if hard != resource.RLIM_INFINITY else min_soft
    if soft < target:
        resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))


def run_snakemake_target(
    target: str | list[str],
    config: Optional[dict] = None,
    configfile: str = "config.yaml",
    dry_run: bool = False,
    cores: int = 4,
    timeout: Optional[int] = None,
) -> dict:
    """
    Replaces pipeline.sh workflow
    """
    targets = [target] if isinstance(target, str) else list(target)
    cmd = ["snakemake"]
    if dry_run:
        cmd.append("-n")
    cmd += targets
    cmd += ["--configfile", configfile, "--cores", str(cores)]
    if config:
        cmd.append("--config")
        cmd += [f"{k}={v}" for k, v in config.items()]

    _ensure_file_descriptor_limit()
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    target_status = {t: Path(t).exists() for t in targets}
    all_exist = all(target_status.values())
    warning_markers = ("Exception", "WorkflowError", "MissingInputException")
    stderr_flagged = any(m in result.stderr for m in warning_markers)
    return {
        "targets": targets,
        "dry_run": dry_run,
        "returncode": result.returncode,
        "succeeded": result.returncode == 0 and (dry_run or all_exist),
        "target_exists_on_disk": target_status,
        "stderr_flagged": stderr_flagged,
        "stdout_tail": result.stdout[-2000:],
        "stderr_tail": result.stderr[-2000:],
    }


def run_deseq2(
    included_samples: list[str],
    sample_attempts: dict[str, int],
    reference_level: str,
    configfile: str = "config.yaml",
    dry_run: bool = False,
) -> dict:
    """
    Request deseq2's outputs via Snakemake
    """
    if not reference_level:
        raise ValueError(
            "reference_level is required. Ground this decision in "
            "get_metadata_conditions() output before calling run_deseq2 -- "
            "do not default to the alphabetically-first Condition value."
        )

    with open(configfile) as fh:
        cfg = yaml.safe_load(fh)

    try:
        resolved_samples = [
            _resolve_run_accession(s, cfg["metadata"]) for s in included_samples
        ]
        resolved_attempts = {
            _resolve_run_accession(s, cfg["metadata"]): a
            for s, a in sample_attempts.items()
        }
    except ValueError as e:
        return {"error": str(e)}

    missing_attempts = set(resolved_samples) - set(resolved_attempts)
    if missing_attempts:
        raise ValueError(
            f"No accepted attempt recorded for: {sorted(missing_attempts)}. "
            "Every included sample needs a QC decision before deseq2 can run."
        )

    config = {
        "included_samples": ",".join(resolved_samples),
        "reference_level": reference_level,
    }
    config.update({f"{s}_attempt": a for s, a in resolved_attempts.items()})

    result = run_snakemake_target(
        ["results/deseq2/deseq2_lfc.txt", "results/deseq2/MA_plot.pdf"],
        config=config,
        configfile=configfile,
        dry_run=dry_run,
    )
    result["reference_level_used"] = reference_level
    return result


def get_top_deseq2_genes(
    results_file: str = "results/deseq2/deseq2_lfc.txt", n: int = 10
) -> dict:
    """
    Reads deseq2's results table and returns the top n rows by adjusted
    p-value as structured data, so the agent can identify genes to
    annotate without guessing at file contents or assuming a standard
    DESeq2 output layout. Row names (gene IDs) come from R's write.table
    convention -- index_col=0 handles the header-row offset that produces.
    """
    df = pd.read_csv(results_file, sep="\t", index_col=0)
    top = df.head(n)
    return {
        "results_file": results_file,
        "n_requested": n,
        "n_returned": len(top),
        "genes": [{"gene_id": gid, **row.to_dict()} for gid, row in top.iterrows()],
    }


def decide_qc_action(qc_summary: dict) -> str:
    """
    Given parse_fastqc_summary()'s output, return "proceed", "retry", or
    "exclude". Each threshold is a specific call, not a guess:
      - No sequences at all: nothing to retry. -> exclude.
      - Per-base quality FAIL: the trim didn't clean the read enough to
        trust alignment. -> retry (with a stricter trimq).
      - Adapter Content FAIL after trimming: same idea,
        trimming under-performed. -> retry.
      - Anything else FAIL/WARN (duplication, GC skew, etc.) without a
        quality or adapter problem: usually reflects real biology or
        library complexity, not something a retrim fixes. -> proceed.
    """
    _MIN_USABLE_DEPTH = 1_000_000
    if qc_summary.get("total_sequences", 0) < _MIN_USABLE_DEPTH:
        return "exclude"
    if qc_summary.get("per_base_quality_status") == "FAIL":
        return "retry"
    if qc_summary.get("adapter_content_status") == "FAIL":
        return "retry"
    return "proceed"


_ANNOTATION_CACHE = {}


def _load_annotation_index(gff_path: str) -> dict:
    if gff_path in _ANNOTATION_CACHE:
        return _ANNOTATION_CACHE[gff_path]

    index = {}
    alias_owners = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = {}
            for kv in fields[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = urllib.parse.unquote(v)

            record = {
                "id": attrs.get("ID"),
                "name": attrs.get("Name"),
                "description": attrs.get("Note", "(no description in CGD annotation)"),
                "classification": attrs.get("orf_classification"),
            }
            gene_key = record["id"] or record["name"]
            aliases = {attrs.get("ID"), attrs.get("Name")}
            aliases.update(a for a in attrs.get("Alias", "").split(",") if a)
            for alias in aliases:
                if not alias:
                    continue
                alias_owners.setdefault(alias, set()).add(gene_key)
                index[alias] = record  # still last-writer-wins for the record itself

    collisions = {a for a, owners in alias_owners.items() if len(owners) > 1}
    _ANNOTATION_CACHE[gff_path] = (index, collisions)
    return index, collisions


def get_gene_annotation(
    gene_id: str,
    gff_path: str = "pipeline/annotation/C_albicans_SC5314_A21_current_features.gff",
) -> dict:
    """
    Exact-match lookup by systematic name or any legacy orf19/CaO19 alias.
    """
    index, collisions = _load_annotation_index(gff_path)
    record = index.get(gene_id)
    if record is None:
        return {
            "found": False,
            "gene_id": gene_id,
            "note": f"'{gene_id}' not present in {gff_path}",
        }
    result = {"found": True, "gene_id": gene_id, **record}
    if gene_id in collisions:
        result["ambiguous"] = True
        result["note"] = (
            f"'{gene_id}' is claimed by multiple genes in {gff_path}; the "
            "returned record is one of several possible matches, not a "
            "unique resolution."
        )
    return result


TOOLS = [
    {
        "name": "run_pipeline_stage",
        "description": (
            "Runs one named stage of the Tag-Seq pipeline: trim, fastqc, "
            "promote_trimmed, align_count, collect_counts, filter_metadata, or "
            "deseq2. Builds the correct Snakemake target path(s) and --config "
            "values automatically -- never construct a raw Snakemake target "
            "path yourself, always call this by stage name instead."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "stage": {
                    "type": "string",
                    "enum": [
                        "trim",
                        "fastqc",
                        "promote_trimmed",
                        "align_count",
                        "collect_counts",
                        "filter_metadata",
                        "deseq2",
                    ],
                },
                "sample": {
                    "type": "string",
                    "description": "Required for trim/fastqc/promote_trimmed/align_count.",
                },
                "attempt": {
                    "type": "integer",
                    "description": "Required for trim/fastqc/promote_trimmed/align_count.",
                },
                "included_samples": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Required for collect_counts.",
                },
                "reference_level": {
                    "type": "string",
                    "description": "Required for deseq2.",
                },
                "dry_run": {"type": "boolean", "default": False},
            },
            "required": ["stage"],
        },
    },
    {
        "name": "run_deseq2",
        "description": "Run differential expression via Snakemake once QC decisions and reference_level are finalized for every included sample.",
        "input_schema": {
            "type": "object",
            "properties": {
                "included_samples": {"type": "array", "items": {"type": "string"}},
                "sample_attempts": {
                    "type": "object",
                    "description": "Map of sample_id -> accepted attempt number. Must cover every entry in included_samples.",
                },
                "reference_level": {
                    "type": "string",
                    "description": "Condition value to use as baseline, e.g. 'WT'. Must be chosen deliberately from get_metadata_conditions() output.",
                },
            },
            "required": ["included_samples", "sample_attempts", "reference_level"],
        },
    },
    {
        "name": "get_top_deseq2_genes",
        "description": (
            "Reads deseq2's results table and returns the top N genes by "
            "adjusted p-value (padj) as structured data. Always use this to "
            "inspect deseq2 output -- never guess at file contents or assume "
            "a standard DESeq2 output format."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "results_file": {
                    "type": "string",
                    "default": "results/deseq2/deseq2_lfc.txt",
                },
                "n": {"type": "integer", "default": 10},
            },
            "required": [],
        },
    },
    {
        "name": "parse_fastqc_summary",
        "description": "Extract PASS/WARN/FAIL module statuses and basic stats from a FastQC report.",
        "input_schema": {
            "type": "object",
            "properties": {
                "qc_dir": {"type": "string"},
                "trimmed_fastq": {"type": "string"},
            },
            "required": ["qc_dir", "trimmed_fastq"],
        },
    },
    {
        "name": "get_metadata_conditions",
        "description": "Read the distinct Condition values and sample IDs from the sample metadata Excel file.",
        "input_schema": {
            "type": "object",
            "properties": {"metadata_xlsx": {"type": "string"}},
            "required": ["metadata_xlsx"],
        },
    },
    {
        "name": "resolve_sample_id",
        "description": (
            "Resolves an SRA run accession or FASTQ filename (e.g. 'SRR18750287') "
            "to its true Sample_ID and Condition from metadata.xlsx. Always call "
            "this before reasoning about which condition or sample an accession "
            "belongs to -- never guess or pattern-match the correspondence."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "fastq_or_srr": {"type": "string", "description": "e.g. 'SRR18750287'"},
                "metadata_xlsx": {
                    "type": "string",
                    "description": "e.g. 'data/metadata.xlsx'",
                },
            },
            "required": ["fastq_or_srr", "metadata_xlsx"],
        },
    },
    {
        "name": "decide_qc_action",
        "description": "Get a recommended QC verdict (proceed/retry/exclude) for a sample from its parsed FastQC summary. You may override this recommendation, but must state why if you do.",
        "input_schema": {
            "type": "object",
            "properties": {
                "qc_summary": {
                    "type": "object",
                    "description": "Output of parse_fastqc_summary()",
                }
            },
            "required": ["qc_summary"],
        },
    },
    {
        "name": "get_gene_annotation",
        "description": "Exact-match lookup of a gene's CGD annotation (name, description, ORF classification) by systematic name or any legacy orf19/CaO19 alias. Returns found=false if not indexed -- treat that as 'no grounded data available', never invent a description.",
        "input_schema": {
            "type": "object",
            "properties": {"gene_id": {"type": "string"}},
            "required": ["gene_id"],
        },
    },
    {
        "name": "verify_claim",
        "description": "Check whether a claimed gene function is actually supported by its retrieved annotation. Call this for every gene claim before including it in your final answer. You may override a failed verdict, but only with an explicit reason grounded in the actual annotation text -- never silently.",
        "input_schema": {
            "type": "object",
            "properties": {
                "gene_id": {"type": "string"},
                "claimed_function": {"type": "string", "description": "Your own-words summary of the gene's function, not a verbatim copy."},
                "annotation_record": {"type": "object", "description": "Full output of get_gene_annotation() for this gene_id."},
            },
            "required": ["gene_id", "claimed_function", "annotation_record"],
        },
    },
]


if __name__ == "__main__":
    print(f"{len(TOOLS)} tool schemas defined:")
    print(json.dumps([t["name"] for t in TOOLS], indent=2))
