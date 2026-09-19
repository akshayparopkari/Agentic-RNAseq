#!/usr/bin/env python3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from evals.cases import ALL_CASES
from evals.harness import Status, run_all, summarize

_STATUS_SYMBOL = {
    Status.PASS: "PASS   ",
    Status.FAIL: "FAIL   ",
    Status.BLOCKED: "BLOCKED",
    Status.ERROR: "ERROR  ",
}


def main() -> int:
    results = run_all(ALL_CASES)

    by_category: dict[str, list] = {}
    for r in results:
        by_category.setdefault(r.category, []).append(r)

    for category, cat_results in by_category.items():
        print(f"\n=== {category} ===")
        for r in cat_results:
            print(f"  [{_STATUS_SYMBOL[r.status]}] {r.case_id:45s} {r.detail}")

    counts = summarize(results)
    total = len(results)
    print(
        f"\n--- {total} cases: "
        f"{counts['PASS']} pass, {counts['FAIL']} fail, "
        f"{counts['BLOCKED']} blocked, {counts['ERROR']} error ---"
    )

    if counts["FAIL"] or counts["ERROR"]:
        print(
            "\nFAIL/ERROR present -- see evals/cases/<id>/output/result.json "
            "for each case's full detail, or evals/README.md for the "
            "documented findings."
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
