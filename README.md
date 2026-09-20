<p align="left">
  <img src="assets/logo.png" alt="Agentic-RNAseq logo" width="200">
</p>

<h1 align="left">Agentic-RNAseq</h1>

**Adaptive QC**: deciding proceed, retry, or exclude from FastQC output.
![Adaptive QC Demo](assets/demo_qc_decision.gif)

**Critic catching an overclaim**: same gene, same annotation, one claim keeps the caveat/hedge (putative) and passes, the other drops it and gets flagged.
![Critic catching an overclaim demo](assets/demo_hedge_drop.gif)

**Grounded reference level**: blocks a silent alphabetical default, forces a deliberate choice from real metadata.
![Grounded reference level demo](assets/demo_reference_level.gif)

This is a weekend project built on top of my published 3'-TagSeq RNA-seq pipeline (github.com/akshayparopkari/RNAseq). The question I wanted to answer: where does an AI agent actually add value on top of a bioinformatics pipeline that already works, without replacing the deterministic tools that make the results reproducible in the first place.

The pipeline itself does not change. The core pipeline workflow of trimming, QC, alignment, counting, and differential expression is exactly as published. What sits on top is an agent that makes the judgment calls a person used to make by hand: whether a sample's QC output is clean enough to align or needs a retrim, which condition should be the reference level in the DESeq2 design, and whether a claim in a generated results summary is actually backed by the gene annotation data from user supplied GFF file.

## Original publication

A Computational Workflow for Analysis of 3' Tag-Seq RNA-seq data.
DOI: https://doi.org/10.1002/cpz1.664

## How it's put together

The agent drives Snakemake one stage at a time instead of kicking off the whole DAG unattended, so it can look at FastQC output before deciding whether alignment should even happen. That QC decision follows fixed rules: a per-base quality or adapter failure means retry with stricter trimming, depth below a usable floor means exclude the sample outright, anything else passes through to alignment.

The part I spent the most time on is the verifier that sits between the DESeq2 results and the written summary. Any claim about what a gene does gets checked against the annotation data actually retrieved for it, not just whether the AI wrote something that sounds right. That check runs in three passes: a cheap similarity filter rules out claims that aren't even about the same topic as the retrieved text, an entailment model checks whether the specific claim actually follows from that text, and a separate rule looks for hedging language in the source that got quietly dropped in the claim, since a model can turn "putative adhesin" into "adhesin" without ever technically contradicting the source.

Alongside all of this is a small eval suite: thirteen cases covering the verifier, the QC decision logic, and how the pipeline wrappers handle bad input, each one stating up front what a correct answer looks like before anything gets run.

## Running it

Needs `snakemake`, `bbduk.sh`, `fastqc`, `STAR`, and `R` with DESeq2 on PATH, `sentence-transformers` installed, and `ANTHROPIC_API_KEY` set.

### Full pipeline, via Snakemake directly

```bash
snakemake <results_dir>/deseq2/deseq2_lfc.txt <results_dir>/deseq2/MA_plot.pdf \
  --configfile config.yaml --cores <cores> \
  --config included_samples=<sample1>,<sample2>,<...> \
  reference_level=<condition_value> \
  <sample1>_attempt=1 <sample2>_attempt=1 <...>
```

Snakemake resolves everything this depends on, trimming, QC, alignment, on its own.

<details>
<summary>Running individual stages instead</summary>

Every target path is built from `results_dir` in config.yaml, `<results_dir>` below stands for whatever you've set that to. Targets are also keyed by the sample's FASTQ-derived accession, not whatever Sample_ID your metadata file uses for that same sample, that resolution only happens inside the Python wrapper, not in Snakemake itself.

One stage, one sample, at a time, useful when you want to look at each stage's output before deciding what to do next, this is exactly the judgment call the agent automates below:

```bash
# trim a sample
snakemake <results_dir>/trim/<sample>_attempt<N>_trimmed.fastq --configfile config.yaml --cores <cores>

# fastqc on that trimmed file
snakemake <results_dir>/qc/<sample>_attempt<N>_trimmed_fastqc.zip --configfile config.yaml --cores <cores>

# once you've reviewed the FastQC output and accepted an attempt, align + count.
# --config <sample>_attempt=<N> is required, nothing guesses which attempt you meant
snakemake <results_dir>/star/<sample>/<sample>ReadsPerGene.out.tab \
  --configfile config.yaml --cores <cores> --config <sample>_attempt=<N>

# once every included sample has an accepted attempt, run DESeq2.
# reference_level is required too, there's no alphabetical default
snakemake <results_dir>/deseq2/deseq2_lfc.txt <results_dir>/deseq2/MA_plot.pdf \
  --configfile config.yaml --cores <cores> \
  --config included_samples=<sample1>,<sample2>,<...> \
  reference_level=<one of your metadata Condition values> \
  <sample1>_attempt=<N1> <sample2>_attempt=<N2> <...>
```
</details>

### Full pipeline, with the agent making the QC and reference-level calls

```python
from agent.planner import run_agent_loop

result = run_agent_loop(
    task=(
        "Process every sample in the metadata through trim, QC, and "
        "alignment, decide proceed/retry/exclude for each based on the "
        "real QC output, choose a reference_level from the real metadata "
        "conditions, run deseq2 once every sample has a final QC decision, "
        "and report the top 10 genes by padj with verified functional "
        "summaries."
    ),
    max_turns=60,
    max_tokens=1024,  # Recommended value. Update to use more or less output tokens
    verbose=True,
)
print(result)
```

Run that from `python`, `ipython`, or a script.

### Tests

```bash
python evals/run_evals.py
```

## Results

Four samples (SRR18750287, SRR18750288, SRR18750289, SRR18750290), reference level WT. 5,651 genes tested, 3 significant at padj < 0.05.

All thirteen eval cases pass. Six cover the verifier, including the two adversarial cases that mattered most to me: a claim that drops a hedge word the annotation used, and a claim that invents a specific mechanism for a real, correctly-identified gene. Three cover the QC decision (clean sample proceeds, adapter failure retries, low depth excludes). Four cover how the pipeline handles bad input: metadata missing a required column, metadata with no rows in it, an invalid pipeline stage name, a missing required argument.

## What worked and what didn't

My first pass at the verifier used plain TF-IDF and cosine similarity, and it topped out around 55% similarity on claims that were nearly word for word identical to the source annotation. Stripping punctuation mangled real biological terms, plural forms were treated as unrelated to their singulars, and filler language diluted the score. Switching to sentence embeddings fixed all of that without any hand-written text cleanup.

But similarity alone turned out not to be enough. A fabricated claim about a real gene stayed close in similarity to the real annotation simply because both were about the same general topic, similarity measures relatedness, not whether one thing actually follows from another. That's what pushed me to add an entailment model as a second pass rather than relying on similarity by itself.

The eval suite earned its keep almost immediately. It caught an inverted condition in the QC depth check that was excluding every clean sample, and a metadata-handling function that silently returned an empty-looking success instead of flagging that a file had no data in it. Both were fixed only because a test written to describe the correct behavior first caught the code doing something else.

One case I deliberately left out for now: a truncated FASTQ file run through the actual trimming tool. Testing that for real needs a genuinely corrupted file and a full run through the pipeline, not a mocked failure, and I'd rather leave it out than fake it.

## What's next

- Building the corrupted/truncated FASTQ test case
- Testing the similarity and entailment thresholds against a wider set of examples than the handful I hand-picked.
- Making model loading reproducible offline instead of depending on a live download every run.
- The entailment check can pass a claim that adds real domain knowledge the annotation never states, it scored "involved in gluconeogenesis" as entailed against text that only says "Phosphoenolpyruvate carboxykinase." True, but exactly the kind of outside knowledge the agent isn't supposed to add. Needs an eval case and a stricter check than entailment alone.
- Getting the agent to verify every retrieved gene, not just the ones with a named function, is inconsistent under prompting alone, a gene with a real function was verified correctly in one run and silently skipped in an identical rerun. Worth enforcing in code rather than depending on the model to remember.
