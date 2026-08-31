"""
Snakefile -- deterministic pipeline stages for the agentic 3' Tag-Seq workflow.

The trim/QC/align commands below are copied verbatim from pipeline.sh in
github.com/akshayparopkari/RNAseq (Paropkari et al., 2023, Current Protocols,
doi: 10.1002/cpz1.664). See ATTRIBUTION.md for what's adapted vs. vendored.

Design notes:
  - `trim` and `fastqc` carry an `attempt` number so every QC retry stays on
    disk as its own file; nothing is silently overwritten.
  - `promote_trimmed` copies whichever attempt was accepted (communicated via
    `--config <sample>_attempt=N`) to a canonical, attempt-free filename.
    Everything downstream depends only on that canonical file.
  - `align_count` does NOT depend on `fastqc`'s output. The QC gate is
    enforced by the agent choosing when to request `promote_trimmed` /
    `align_count`, not by Snakemake's dependency graph.
  - `collect_counts` only aggregates whichever samples are explicitly listed
    in `--config included_samples=a,b,c`, so an excluded sample is simply
    never in that list.
  - `filter_metadata` derives inclusion purely from gene_raw_counts.txt's own
    column headers -- one source of truth, no separate state file needed.
  - `deseq2` requires an explicit `reference_level` config value. It does
    NOT fall through to deseq.R's own alphabetical-default behavior,
    since that default is not guaranteed to be the biologically correct
    baseline -- see the note in deseq.R.
"""

import glob

import pandas as pd

configfile: "config.yaml"

wildcard_constraints:
    attempt=r"\d+"


def get_samples():
    """Sample IDs, read once from the metadata file's Sample_ID column."""
    meta = pd.read_excel(config["metadata"])
    return sorted(meta["Sample_ID"].dropna().astype(str).tolist())


SAMPLES = get_samples()


def find_raw_fastq(wildcards):
    """
    Locate the one raw FASTQ file for a sample. Matches either an exact
    accession filename (SRR18750287.fastq.gz) or the original pipeline's
    sampleid_restofname.fastq(.gz) convention.
    """
    matches = (
        glob.glob(f"{config['raw_fastq_dir']}/{wildcards.sample}.fastq*")
        + glob.glob(f"{config['raw_fastq_dir']}/{wildcards.sample}_*.fastq*")
    )
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one raw FASTQ file for sample '{wildcards.sample}' "
            f"in {config['raw_fastq_dir']}, found {len(matches)}: {matches}"
        )
    return matches[0]


def accepted_attempt_path(wildcards):
    """
    Which attempt was accepted for this sample, read from a required
    per-sample config override, e.g. --config sample1_attempt=2
    No silent default -- an unset value is a real error.
    """
    key = f"{wildcards.sample}_attempt"
    if key not in config:
        raise ValueError(
            f"No accepted attempt recorded for sample '{wildcards.sample}'. "
            f"Pass --config {key}=<N> once the QC decision has been made."
        )
    attempt = config[key]
    return f"{config['results_dir']}/trim/{wildcards.sample}_attempt{attempt}_trimmed.fastq"


def required_reference_level(wildcards):
    if "reference_level" not in config:
        raise ValueError(
            "No reference_level set. Pass --config reference_level=<value> "
            "explicitly -- this must be a deliberate choice, not a default."
        )
    return config["reference_level"]


def trim_quality_threshold(wildcards):
    """
    Per-(sample, attempt) trimq override, e.g. --config SRR123_attempt2_trimq=30
    for a stricter retry after a failed QC check. Defaults to 20, the
    original published pipeline's value, when no override is given.
    """
    key = f"{wildcards.sample}_attempt{wildcards.attempt}_trimq"
    return config.get(key, 20)


rule all:
    # Convenience target for a manual, full run. The agent never requests
    # this directly -- it requests the intermediate targets below, one
    # decision at a time.
    input:
        f"{config['results_dir']}/deseq2/deseq2_lfc.txt",
        f"{config['results_dir']}/deseq2/MA_plot.pdf",


rule trim:
    """Adapter + polyA trimming. Mirrors CMD1 in pipeline.sh exactly."""
    input:
        find_raw_fastq
    output:
        f"{config['results_dir']}/trim/{{sample}}_attempt{{attempt}}_trimmed.fastq"
    params:
        truseq_ref=config["truseq_ref"],
        polya_ref=config["polya_ref"],
        trimq=trim_quality_threshold
    log:
        f"{config['results_dir']}/logs/trim_{{sample}}_attempt{{attempt}}.log"
    shell:
        "bbduk.sh in={input} out={output} "
        "ref={params.truseq_ref},{params.polya_ref} "
        "k=13 ktrim=r mink=5 qtrim=r trimq={params.trimq} minlength=20 "
        "> {log} 2>&1"


rule fastqc:
    """Read quality report on a trimmed FASTQ. Mirrors CMD2 exactly."""
    input:
        f"{config['results_dir']}/trim/{{sample}}_attempt{{attempt}}_trimmed.fastq"
    output:
        f"{config['results_dir']}/qc/{{sample}}_attempt{{attempt}}_trimmed_fastqc.zip"
    params:
        outdir=f"{config['results_dir']}/qc",
        threads=config["threads"]
    log:
        f"{config['results_dir']}/logs/fastqc_{{sample}}_attempt{{attempt}}.log"
    shell:
        "fastqc -t {params.threads} --nogroup -o {params.outdir} {input} > {log} 2>&1"


rule promote_trimmed:
    """
    Copies the accepted attempt's trimmed FASTQ to a canonical, attempt-free
    filename. Everything downstream depends only on this file, never on a
    specific attempt number directly.
    """
    input:
        accepted_attempt_path
    output:
        f"{config['results_dir']}/trim/canonical/{{sample}}_trimmed.fastq"
    log:
        f"{config['results_dir']}/logs/promote_{{sample}}.log"
    shell:
        "cp {input} {output} > {log} 2>&1"


rule align_count:
    """
    Alignment + gene counting. Mirrors CMD3 in pipeline.sh exactly. The
    output prefix MUST be exactly the sample ID with nothing appended --
    format_counts_table.py extracts sample IDs with a naive filename
    string-replace and has no tolerance for anything extra.
    """
    input:
        f"{config['results_dir']}/trim/canonical/{{sample}}_trimmed.fastq"
    output:
        counts=f"{config['results_dir']}/star/{{sample}}/{{sample}}ReadsPerGene.out.tab",
        bam=f"{config['results_dir']}/star/{{sample}}/{{sample}}Aligned.sortedByCoord.out.bam"
    params:
        prefix=f"{config['results_dir']}/star/{{sample}}/{{sample}}",
        threads=config["threads"],
        genome_dir=config["genome_dir"]
    log:
        f"{config['results_dir']}/logs/align_{{sample}}.log"
    shell:
        "STAR --runThreadN {params.threads} --genomeDir {params.genome_dir} "
        "--readFilesIn {input} --outFilterType BySJout --outFilterMultimapNmax 25 "
        "--alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverLmax 0.3 --alignIntronMin 20 --alignIntronMax 1000000 "
        "--alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD "
        "--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts "
        "--outFileNamePrefix {params.prefix} > {log} 2>&1"


rule collect_counts:
    """
    Aggregates every requested sample's counts into one matrix via the
    original format_counts_table.py, unmodified. Which samples are
    aggregated is fully determined by --config included_samples=a,b,c --
    an excluded sample is simply never in that list.
    """
    input:
        lambda wc: expand(
            f"{config['results_dir']}/star/{{sample}}/{{sample}}ReadsPerGene.out.tab",
            sample=config["included_samples"].split(",")
        )
    output:
        f"{config['results_dir']}/counts/gene_raw_counts.txt"
    params:
        star_dir=f"{config['results_dir']}/star",
        out_dir=f"{config['results_dir']}/counts"
    log:
        f"{config['results_dir']}/logs/collect_counts.log"
    shell:
        "python3 pipeline/format_counts_table.py {params.star_dir} -o {params.out_dir} > {log} 2>&1"


rule filter_metadata:
    """
    Keeps only metadata rows for samples that actually made it into
    gene_raw_counts.txt. No separate inclusion list -- the count matrix's
    own columns are the single source of truth.
    """
    input:
        metadata=config["metadata"],
        counts=f"{config['results_dir']}/counts/gene_raw_counts.txt"
    output:
        f"{config['results_dir']}/metadata_filtered.xlsx"
    log:
        f"{config['results_dir']}/logs/filter_metadata.log"
    script:
        "scripts/filter_metadata.py"


rule deseq2:
    """Differential expression. Calls deseq.R -- see ATTRIBUTION.md."""
    input:
        counts=f"{config['results_dir']}/counts/gene_raw_counts.txt",
        metadata=f"{config['results_dir']}/metadata_filtered.xlsx"
    output:
        results=f"{config['results_dir']}/deseq2/deseq2_lfc.txt",
        ma_plot=f"{config['results_dir']}/deseq2/MA_plot.pdf"
    params:
        reference_level=required_reference_level
    log:
        f"{config['results_dir']}/logs/deseq2.log"
    shell:
        "Rscript --vanilla pipeline/deseq.R {input.counts} {input.metadata} "
        "{output.results} {output.ma_plot} {params.reference_level} > {log} 2>&1"
