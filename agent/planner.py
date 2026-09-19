import json

from anthropic import Anthropic

from agent.tools import (
    TOOLS,
    run_deseq2,
    decide_qc_action,
    resolve_sample_id,
    run_pipeline_stage,
    get_gene_annotation,
    get_top_deseq2_genes,
    parse_fastqc_summary,
    get_metadata_conditions,
)
from agent.critic import verify_claim_against_annotation

MODEL = "claude-haiku-4-5-20251001"

SYSTEM_PROMPT = """You are orchestrating a 3' Tag-Seq RNA-seq pipeline via Snakemake tool calls.

Every task message defines the complete scope of what you are authorized to
do this turn. The numbered steps below describe the overall pipeline's logic
and order so you understand how the stages relate to each other -- they are
not standing authorization to act on every sample, stage, or pipeline step
whenever you're invoked. Only act on the sample(s) or stage(s) the current
task actually names. Resolving what a task refers to is always in scope,
even when it doesn't spell out literal IDs: if a task says "all four
samples" or names a run accession without its Sample_ID, call the
appropriate lookup tool (get_metadata_conditions, resolve_sample_id) to
ground it before acting -- do not ask the user for information a tool can
already give you. Only ask the user, or stop and report, when no available
tool can resolve what's being referred to. If you notice other samples are
incomplete or pending while working, mention it in your final summary, but
do not call tools for anything the task didn't explicitly ask for. Once
you've completed everything the current task asked for, stop and report --
do not look for additional work to perform.

Your job, in order:
1. Call get_metadata_conditions to see the real Condition values and sample IDs.
   Choose reference_level deliberately from what's actually there -- never
   assume "WT" or alphabetical order. State your reasoning before choosing.
2. For each sample named in your current task: request its trim + fastqc
   targets via run_pipeline_stage, call parse_fastqc_summary on the result,
   then call decide_qc_action for a recommended verdict. You may override
   its recommendation, but only with an explicit reason grounded in the
   actual QC numbers -- never silently.
   - "proceed": request that sample's align_count target.
   - "retry": request another attempt with a stricter trimq (pass
     {sample}_attempt{N}_trimq in config), then re-check QC. Track retries
     yourself -- after 2 with no improvement, exclude instead of looping.
   - "exclude": skip align_count for this sample. Note why.
3. Once every included sample has a final result or is excluded, and only
   when the current task explicitly asks you to, call run_deseq2 with the
   included sample list, their accepted attempts, and your chosen
   reference_level.
4. For each gene you report on, call get_gene_annotation(gene_id). If
   found=false, say plainly that no annotation was found for that gene --
   do not describe its function from your own knowledge, and do not call
   verify_claim for it. If found=true, a real record was retrieved --
   never describe that gene as having no annotation, not annotated, or
   unavailable, no matter how sparse the annotation is (even "protein of
   unknown function" is a real, retrieved record, not an absence of one)
   and no matter what you decide to do with verify_claim for it. Summarize
   its function in your own words -- do not copy the description verbatim,
   a verbatim copy defeats the point of this step. This applies even when
   the annotation is sparse or when your summary only adds a detail like
   where or when the gene is induced -- call verify_claim for every
   found=true gene without exception. There is no claim too small or too
   well-supported-sounding to skip this step for. Call verify_claim with
   that gene_id, your own-words summary as claimed_function, and the full
   get_gene_annotation result as annotation_record. Only include the claim
   in your final report if verify_claim returns supported=true, or if you
   override a failed verdict with an explicit reason grounded in the
   actual annotation text -- never silently skip verify_claim, and never
   silently keep a claim it failed. Summarize statistics (fold change,
   p-value) directly from run_deseq2's output either way. Report each gene
   as one compact line: gene_id, log2FoldChange, padj, then your verified
   summary of its function. No headers, no per-gene bolded titles, no
   restating the question. Facts only.

Your final answer must restate every decision made this session, not just
the most recent step. Never infer, pattern-match, or assume a correspondence
between two identifiers referring to the same entity (sample, gene, file,
or otherwise). Any such correspondence must come from an explicit tool result.
The same rule applies to file paths and Snakemake targets: never construct a
raw target path yourself, always call run_pipeline_stage with the stage name
and its required arguments (sample, attempt, included_samples, or
reference_level) and let the tool derive the real path from config.yaml and
the Snakefile's own conventions. If no tool can resolve something, say so
plainly and stop rather than proceeding on an assumption -- this includes
stopping rather than guessing at an identifier (e.g. inferring a nearby
accession number) when a lookup tool fails to resolve one you need. Every
tool result is ground truth. If a tool call fails, stop and report what
failed rather than guessing what would have happened.
"""

TOOL_IMPL = {
    "run_pipeline_stage": run_pipeline_stage,
    "run_deseq2": run_deseq2,
    "parse_fastqc_summary": parse_fastqc_summary,
    "get_metadata_conditions": get_metadata_conditions,
    "decide_qc_action": decide_qc_action,
    "get_gene_annotation": get_gene_annotation,
    "resolve_sample_id": resolve_sample_id,
    "get_top_deseq2_genes": get_top_deseq2_genes,
    "verify_claim": lambda gene_id, claimed_function, annotation_record: (
        verify_claim_against_annotation(
            gene_id, claimed_function, annotation_record
        ).as_dict()
    ),
}


def run_agent_loop(
    task: str, max_turns: int = 20, max_tokens: int = 512, verbose: bool = False
) -> str:
    """
    Drives one Claude conversation to completion, dispatching tool calls to
    the real functions in agent/tools.py. Returns the full text the agent
    produced across every turn, not just the final one, so intermediate
    determinations (e.g. reference_level) aren't dropped from the result.
    """
    client = Anthropic()
    messages = [{"role": "user", "content": task}]
    transcript = []  # every non-empty text block the agent produced, in order

    for turn in range(max_turns):
        response = client.messages.create(
            model=MODEL,
            max_tokens=max_tokens,
            system=SYSTEM_PROMPT,
            tools=TOOLS,
            messages=messages,
        )
        messages.append({"role": "assistant", "content": response.content})

        for block in response.content:
            if block.type == "text" and block.text.strip():
                transcript.append(block.text.strip())
                if verbose:
                    print(f"[turn {turn}] TEXT: {block.text}")
            elif block.type == "tool_use" and verbose:
                print(f"[turn {turn}] CALL: {block.name}({block.input})")

        if response.stop_reason != "tool_use":
            return "\n\n".join(transcript)

        if response.stop_reason == "max_tokens":
            if verbose:
                print(f"[turn {turn}] WARNING: response truncated by max_tokens")
            transcript.append(
                "[WARNING: this response was truncated by the max_tokens limit and is incomplete]"
            )
            return "\n\n".join(transcript)
        elif response.stop_reason != "tool_use":
            return "\n\n".join(transcript)

        tool_results = []
        for block in response.content:
            if block.type != "tool_use":
                continue
            func = TOOL_IMPL.get(block.name)
            try:
                output = (
                    func(**block.input)
                    if func
                    else {"error": f"unknown tool {block.name}"}
                )
            except Exception as e:
                output = {
                    "error": f"{type(e).__name__}: {e}"
                }  # report to agent, don't crash the loop
            if verbose:
                print(f"[turn {turn}] RESULT {block.name}: {output}")
            tool_results.append(
                {
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": json.dumps(output, default=str),
                }
            )
        messages.append({"role": "user", "content": tool_results})

    raise RuntimeError(f"Agent loop did not reach a final answer in {max_turns} turns")
