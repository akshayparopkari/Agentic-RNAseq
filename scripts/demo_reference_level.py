# demo_reference_level.py
from pprint import pprint
from agent.tools import get_metadata_conditions, run_deseq2

conditions = get_metadata_conditions("data/metadata.xlsx")
print("\n--- Conditions found in metadata ---")
print(conditions)
print(
    f"\nNote: sorted alphabetically, {conditions['condition_levels'][0]!r} comes "
    f"first. That is NOT automatically used as the reference level."
)

print("\n--- Attempting run_deseq2 without a reference_level ---")
try:
    run_deseq2(
        included_samples=conditions["sample_ids"],
        sample_attempts={s: 1 for s in conditions["sample_ids"]},
        reference_level="",
    )
except ValueError as e:
    print(f"Blocked: {e}")

print("\n--- Grounding the choice deliberately in the real metadata (dry run) ---")
result = run_deseq2(
    included_samples=conditions["sample_ids"],
    sample_attempts={s: 1 for s in conditions["sample_ids"]},
    reference_level="WT",  # the wild-type baseline, chosen deliberately -- not alphabetical
    dry_run=True,
)
print(f"succeeded: {result['succeeded']}, targets: {result['targets']}")
print()