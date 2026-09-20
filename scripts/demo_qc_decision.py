# demo_qc_decision.py
from pprint import pprint
from agent.tools import decide_qc_action

print("\n--- Case 1: clean sample, high depth ---")
clean = {
    "total_sequences": 3_183_204,
    "adapter_content_status": "PASS",
    "per_base_quality_status": "PASS",
    "duplication_status": "WARN",  # benign for 3' Tag-Seq, not a retrim trigger
}
pprint(f"verdict: {decide_qc_action(clean)}")

print("\n--- Case 2: adapter content failed after trim ---")
noisy = {
    "total_sequences": 2_754_310,
    "adapter_content_status": "FAIL",
    "per_base_quality_status": "PASS",
    "duplication_status": "PASS",
}
pprint(f"verdict: {decide_qc_action(noisy)}")

print("\n--- Case 3: depth below the usable floor ---")
low_depth = {
    "total_sequences": 612_040,
    "adapter_content_status": "PASS",
    "per_base_quality_status": "PASS",
    "duplication_status": "PASS",
}
pprint(f"verdict: {decide_qc_action(low_depth)}")
