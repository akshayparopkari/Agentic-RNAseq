# demo_hedge_drop.py
from pprint import pprint
from agent.critic import verify_claim_against_annotation

annotation_record = {
    "found": True,
    "gene_id": "orf19.5741",
    "description": (
        "ALS1: putative GPI-anchored adhesin, inferred from homology to "
        "S. cerevisiae ALS family; implicated in biofilm formation."
    ),
}

result = verify_claim_against_annotation(
    "orf19.5741",
    "ALS1 is a putative adhesin, inferred from homology, implicated in biofilm formation.",
    annotation_record,
)
print("\n--- Claim 1: honest, hedge preserved ---")
pprint(result.as_dict())

print("\n--- Claim 2: same gene, hedge dropped ---")
result = verify_claim_against_annotation(
    "orf19.5741",
    "ALS1 is an adhesin implicated in biofilm formation.",
    annotation_record,
)
pprint(result.as_dict())
