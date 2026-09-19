#!/usr/bin/env python3
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

from agent import tools, critic
from evals.harness import Case

CASES_DIR = Path(__file__).parent / "cases"


# CRITIC/VERIFIER layer case check


def _case_hallucinated_gene_id():
    return critic.verify_claim_against_annotation(
        gene_id="orf19.5741",
        claimed_function="Encodes a cell-wall adhesin implicated in biofilm formation.",
        annotation_record={
            "found": True,
            "gene_id": "orf19.1816",
            "name": "ALS3",
            "description": (
                "cell-wall adhesin, GPI-anchored, implicated in adhesion "
                "and biofilm formation. Paralog of ALS1."
            ),
        },
    )


def _case_empty_annotation_retrieved():
    return critic.verify_claim_against_annotation(
        gene_id="orf19.5741",
        claimed_function="Encodes a cell-wall adhesin implicated in biofilm formation.",
        annotation_record={
            "found": False,
            "gene_id": "orf19.5741",
            "note": "'orf19.5741' not present in the annotation file",
        },
    )


def _case_empty_claim():
    return critic.verify_claim_against_annotation(
        gene_id="orf19.5741",
        claimed_function="   ",
        annotation_record={
            "found": True,
            "gene_id": "orf19.5741",
            "name": "ALS1",
            "description": "cell-wall adhesin, GPI-anchored.",
        },
    )


def _case_grounded_claim_should_pass():
    """
    Positive check for claimed function of a gene. Must pass.
    """
    return critic.verify_claim_against_annotation(
        gene_id="orf19.5741",
        claimed_function="Encodes a cell-wall adhesin.",
        annotation_record={
            "found": True,
            "gene_id": "orf19.5741",
            "name": "ALS1",
            "description": (
                "cell-wall adhesin, GPI-anchored, mediates adhesion to "
                "host epithelial cells."
            ),
        },
    )


def _case_dropped_hedge_overclaim():
    """
    Hedge check, dropping putative/inferred qualifiers in claimed function of a gene
    is misleading and must be reported in the claim.
    """
    return critic.verify_claim_against_annotation(
        gene_id="orf19.2475",
        claimed_function="Is a cell-wall adhesin.",
        annotation_record={
            "found": True,
            "gene_id": "orf19.2475",
            "name": None,
            "description": (
                "putative cell-wall adhesin (inferred from homology to S. "
                "cerevisiae FLO11); function not experimentally validated "
                "in C. albicans."
            ),
        },
    )


def _case_hallucinated_mechanism():
    """
    Is the claimed function imagined/hallucinated by the agent?
    """
    return critic.verify_claim_against_annotation(
        gene_id="orf19.1774",
        claimed_function=(
            "Is transcriptionally repressed by Efg1 under hyphal-inducing "
            "conditions."
        ),
        annotation_record={
            "found": True,
            "gene_id": "orf19.1774",
            "name": "HWP1",
            "description": (
                "hyphal cell wall protein, GPI-anchored, induced during "
                "hyphal morphogenesis; established marker gene for the "
                "yeast-to-hyphae transition."
            ),
        },
    )


CRITIC_CASES = [
    Case(
        id="hallucinated_gene_id",
        category="critic",
        description="Claim's gene_id never appears in the retrieved annotation text.",
        run=_case_hallucinated_gene_id,
        expected=False,
        matches=lambda actual, expected: actual.supported == expected
        and actual.check == "gene_id_retrieved",
    ),
    Case(
        id="empty_annotation_retrieved",
        category="critic",
        description="Annotation retrieval came back empty for this gene_id.",
        run=_case_empty_annotation_retrieved,
        expected=False,
        matches=lambda actual, expected: actual.supported == expected
        and actual.check == "gene_id_retrieved",
    ),
    Case(
        id="empty_claim",
        category="critic",
        description="Reporter step produced a blank/whitespace-only claim.",
        run=_case_empty_claim,
        expected=False,
        matches=lambda actual, expected: actual.supported == expected
        and actual.check == "claim_nonempty",
    ),
    Case(
        id="grounded_claim_should_pass",
        category="critic",
        description="Faithful, accurate restatement of the retrieved annotation.",
        run=_case_grounded_claim_should_pass,
        expected=True,
        matches=lambda actual, expected: actual.supported == expected,
        blocked_on="agent.critic._check_semantic_support",
        notes="Simplest possible semantic case -- good first target once Block 5 is written.",
    ),
    Case(
        id="dropped_hedge_overclaim",
        category="critic",
        description="Claim drops the annotation's 'putative / inferred from homology' hedge and states it as fact.",
        run=_case_dropped_hedge_overclaim,
        expected=False,
        matches=lambda actual, expected: actual.supported == expected,
        blocked_on="agent.critic._check_semantic_support",
        notes="The flagship grounding-gap case -- see critic.py module docstring design question 3.",
    ),
    Case(
        id="hallucinated_mechanism",
        category="critic",
        description="Correct gene, but claim states a specific mechanism the annotation never mentions.",
        run=_case_hallucinated_mechanism,
        expected=False,
        matches=lambda actual, expected: actual.supported == expected,
        blocked_on="agent.critic._check_semantic_support",
    ),
]

_CLEAN_QC_SUMMARY = {
    "modules": {
        "Per base sequence quality": "PASS",
        "Sequence Duplication Levels": "PASS",
        "Overrepresented sequences": "PASS",
    },
    "total_sequences": 8_500_000,
    "pct_gc": 45,
    "sequence_length": "50",
    "adapter_content_status": "PASS",
    "duplication_status": "PASS",
    "per_base_quality_status": "PASS",
}

_NOISY_BUT_RECOVERABLE_QC_SUMMARY_ATTEMPT1 = {
    **_CLEAN_QC_SUMMARY,
    "modules": {**_CLEAN_QC_SUMMARY["modules"], "Overrepresented sequences": "FAIL"},
    "adapter_content_status": "FAIL",
}

_LOW_DEPTH_QC_SUMMARY = {
    **_CLEAN_QC_SUMMARY,
    "total_sequences": 120_000,
}


def _case_clean_sample_should_proceed():
    return tools.decide_qc_action(_CLEAN_QC_SUMMARY)


def _case_noisy_sample_should_retrim():
    return tools.decide_qc_action(_NOISY_BUT_RECOVERABLE_QC_SUMMARY_ATTEMPT1)


def _case_low_depth_should_exclude():
    return tools.decide_qc_action(_LOW_DEPTH_QC_SUMMARY)


QC_DECISION_CASES = [
    Case(
        id="clean_sample_should_proceed",
        category="qc_decision",
        description="All modules PASS, healthy read depth -- should proceed straight to alignment.",
        run=_case_clean_sample_should_proceed,
        expected="proceed",
        blocked_on="agent.tools.decide_qc_action",
    ),
    Case(
        id="noisy_sample_should_retrim",
        category="qc_decision",
        description="Adapter content WARN on attempt 1, everything else clean -- recoverable, should retrim rather than exclude outright.",
        run=_case_noisy_sample_should_retrim,
        expected="retry",
        blocked_on="agent.tools.decide_qc_action",
        notes="Pair with a second call using a clean attempt-2 summary once Block 3 exists, "
        "to demonstrate the retrim->proceed recovery path end to end.",
    ),
    Case(
        id="low_depth_should_exclude",
        category="qc_decision",
        description="Total sequences far below a usable depth floor -- retrimming can't fix low input material, should exclude.",
        run=_case_low_depth_should_exclude,
        expected="exclude",
        blocked_on="agent.tools.decide_qc_action",
    ),
]


def _case_malformed_metadata_missing_condition():
    return tools.get_metadata_conditions(
        str(
            CASES_DIR
            / "malformed_metadata_missing_condition_column"
            / "fixtures"
            / "metadata.xlsx"
        )
    )


def _case_malformed_metadata_empty():
    return tools.get_metadata_conditions(
        str(CASES_DIR / "malformed_metadata_empty_file" / "fixtures" / "metadata.xlsx")
    )


PIPELINE_STAGE_CASES = [
    Case(
        id="malformed_metadata_missing_condition_column",
        category="pipeline_stage",
        description="metadata.xlsx has Sample_ID and Run_Accession but no Condition column.",
        run=_case_malformed_metadata_missing_condition,
        expected={"error": "metadata file has no 'Condition' column"},
        notes="Real fixture file, real pandas read -- no mocking.",
    ),
    Case(
        id="malformed_metadata_empty_file",
        category="pipeline_stage",
        description="metadata.xlsx has correct headers but zero data rows.",
        run=_case_malformed_metadata_empty,
        expected={"error": "metadata file has zero data rows"},
        matches=lambda actual, expected: isinstance(actual, dict) and "error" in actual,
        notes="Currently FAILS -- see evals/README.md 'Findings' for why this is a real gap "
        "worth fixing, not a mistake in the eval case.",
    ),
    Case(
        id="unknown_stage_name",
        category="pipeline_stage",
        description="run_pipeline_stage called with a stage name that isn't valid.",
        run=lambda: tools.run_pipeline_stage(stage="nonexistent"),
        expected={"error": "..."},
        matches=lambda a, e: isinstance(a, dict) and "error" in a,
    ),
    Case(
        id="missing_required_stage_arg",
        category="pipeline_stage",
        description="run_pipeline_stage called for 'trim' without sample/attempt.",
        run=lambda: tools.run_pipeline_stage(stage="trim"),
        expected={"error": "..."},
        matches=lambda a, e: isinstance(a, dict) and "error" in a,
    ),
]


ALL_CASES = CRITIC_CASES + QC_DECISION_CASES + PIPELINE_STAGE_CASES

if __name__ == "__main__":
    print(f"{len(ALL_CASES)} cases defined:")
    for c in ALL_CASES:
        flag = f" [blocked_on={c.blocked_on}]" if c.blocked_on else ""
        print(f"  [{c.category}] {c.id}{flag}")
