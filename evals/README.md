# Eval suite

Sprint plan Block 6. 13 cases across the three places this project makes a
judgment call or wraps an external tool: `agent/critic.py`, the still-unwritten
`agent/tools.decide_qc_action` (Block 3), and the error handling in
`agent/tools.py`'s pipeline-stage wrappers.

Run from the repo root:

```
python evals/run_evals.py
```

## Current results

Last run in this cloud sandbox (no bioinformatics tools installed — bbduk,
FastQC, STAR, Snakemake are all absent here, and nothing in the QC-decision
or semantic-verification code paths has been written yet):

```
13 cases: 4 pass, 1 fail, 7 blocked, 1 error
```

| status | meaning |
|---|---|
| PASS | ran against real code, matched what a correct implementation should produce |
| FAIL | ran against real code, did **not** match — an actual gap worth fixing |
| BLOCKED | didn't run at all, because it depends on something declared in advance as unfinished (Block 3, Block 5's semantic check) or a binary this sandbox doesn't have (bbduk, STAR) |
| ERROR | raised something that was **not** declared in advance — a surprise |

The distinction between BLOCKED and ERROR is deliberate. Every BLOCKED case
sets `blocked_on` up front, before it's ever run, naming exactly what it's
waiting on. An ERROR means a case hit a failure mode nobody anticipated when
it was written — that's the one status here that means "stop and look."

## Findings

Two real gaps this suite surfaced in the existing pipeline-stage wrappers
(`agent/tools.py`), independent of the two blocks still pending:

- **`malformed_metadata_empty_file` — FAIL.** `get_metadata_conditions()`
  already returns a clean `{"error": ...}` dict when the `Condition` column
  is missing entirely, but a metadata file with the right headers and zero
  data rows silently returns `{"condition_levels": [], "sample_count": 0,
  "sample_ids": []}` instead of also being treated as an error. Worth
  hardening: a downstream planning step that doesn't explicitly check
  `sample_count == 0` would happily proceed with no samples.
- **`missing_star_binary` — ERROR.** None of the stage wrappers
  (`trim_sample`, `qc_sample`, `align_and_count`, `run_deseq2`) check
  whether their underlying binary is even on `PATH` before calling
  `subprocess.run`. When it isn't, the call raises an uncaught
  `FileNotFoundError` instead of returning the same structured
  `{"returncode": ..., "stderr_tail": ...}` shape every other failure mode
  in this file uses. An agent loop dispatching tool calls would see a raw
  Python exception here instead of a result it can reason about — worth
  wrapping every `subprocess.run` call in the same try/except and folding
  "binary not found" into the existing returncode-based error contract.

Both are logic/design gaps, not sandbox artifacts — they'll reproduce
identically in an environment that does have STAR and real metadata files.

## What's genuinely untestable in this sandbox

`truncated_fastq_trim_stage` needs bbduk.sh actually installed to observe
real behavior on a corrupted `.fastq.gz` (does it recover partial reads,
fail cleanly, fail with a confusing error?). This sandbox has none of the
bioinformatics binaries installed at all, so this case can't even get past
"is the tool there" to ask the real question. Needs to run in an environment
with the pipeline's actual dependencies (e.g. the `rnaseq` mamba env from
the original working session, or a machine with the tools installed) before
it means anything.

## What flips BLOCKED to PASS/FAIL

- Writing `agent/tools.decide_qc_action` (Block 3) unblocks
  `clean_sample_should_proceed`, `noisy_sample_should_retrim`, and
  `low_depth_should_exclude` — all three call it with realistic
  `parse_fastqc_summary()`-shaped dicts, no real FastQC run required.
- Writing `agent/critic._check_semantic_support` (Block 5's core logic)
  unblocks `grounded_claim_should_pass`, `dropped_hedge_overclaim`, and
  `hallucinated_mechanism`. The three mechanical critic checks
  (`hallucinated_gene_id`, `empty_annotation_retrieved`, `empty_claim`)
  already pass today and aren't affected by Block 5.
- `truncated_fastq_trim_stage` needs a real BBDuk install, not more code.

## Case list

| id | category | tests |
|---|---|---|
| `hallucinated_gene_id` | critic | claim cites a gene_id absent from the retrieved annotation |
| `empty_annotation_retrieved` | critic | annotation lookup came back empty |
| `empty_claim` | critic | blank claimed_function |
| `grounded_claim_should_pass` | critic | faithful, accurate restatement of the annotation |
| `dropped_hedge_overclaim` | critic | claim drops a "putative / inferred from homology" hedge |
| `hallucinated_mechanism` | critic | correct gene, invented mechanism not in the annotation |
| `clean_sample_should_proceed` | qc_decision | all FastQC modules PASS, healthy depth |
| `noisy_sample_should_retrim` | qc_decision | one WARN, otherwise clean — recoverable |
| `low_depth_should_exclude` | qc_decision | total_sequences far below a usable floor |
| `malformed_metadata_missing_condition_column` | pipeline_stage | metadata.xlsx missing the Condition column |
| `malformed_metadata_empty_file` | pipeline_stage | metadata.xlsx with headers but zero rows |
| `missing_star_binary` | pipeline_stage | align_and_count when STAR isn't on PATH |
| `truncated_fastq_trim_stage` | pipeline_stage | trim_sample on a corrupted FASTQ (needs real BBDuk) |
