#!/usr/bin/env python3
import json
import time
import traceback
from enum import Enum
from typing import Any, Callable, Optional
from pathlib import Path
from dataclasses import dataclass

EVALS_DIR = Path(__file__).parent


class Status(str, Enum):
    PASS = "PASS"
    FAIL = "FAIL"
    BLOCKED = "BLOCKED"
    ERROR = "ERROR"


def _default_matches(actual: Any, expected: Any) -> bool:
    return actual == expected


@dataclass
class Case:
    id: str
    category: str
    description: str
    run: Callable[[], Any]
    expected: Any
    matches: Callable[[Any, Any], bool] = _default_matches
    blocked_on: Optional[str] = None
    notes: str = ""


@dataclass
class CaseResult:
    case_id: str
    category: str
    status: Status
    detail: str
    actual: Any = None
    expected: Any = None
    duration_ms: float = 0.0


def _to_jsonable(value: Any) -> Any:
    if hasattr(value, "as_dict"):
        return value.as_dict()
    if isinstance(value, list):
        return [_to_jsonable(v) for v in value]
    return value


def _write_result(out_dir: Path, result: CaseResult, case: Case) -> None:
    payload = {
        "case_id": result.case_id,
        "category": result.category,
        "status": result.status.value,
        "detail": result.detail,
        "actual": _to_jsonable(result.actual),
        "expected": _to_jsonable(result.expected),
        "duration_ms": round(result.duration_ms, 2),
        "description": case.description,
        "blocked_on": case.blocked_on,
        "notes": case.notes,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "result.json").write_text(json.dumps(payload, indent=2, default=str))


def run_case(case: Case) -> CaseResult:
    out_dir = EVALS_DIR / "cases" / case.id / "output"
    start = time.time()

    try:
        actual = case.run()
    except (NotImplementedError, FileNotFoundError) as e:
        duration = (time.time() - start) * 1000
        if case.blocked_on:
            result = CaseResult(
                case_id=case.id,
                category=case.category,
                status=Status.BLOCKED,
                detail=f"Blocked on {case.blocked_on}: {e}",
                expected=case.expected,
                duration_ms=duration,
            )
        else:
            result = CaseResult(
                case_id=case.id,
                category=case.category,
                status=Status.ERROR,
                detail=f"Unanticipated {type(e).__name__} (no blocked_on set "
                f"for this case, so this wasn't expected): {e}",
                expected=case.expected,
                duration_ms=duration,
            )
        _write_result(out_dir, result, case)
        return result
    except Exception:
        duration = (time.time() - start) * 1000
        result = CaseResult(
            case_id=case.id,
            category=case.category,
            status=Status.ERROR,
            detail=f"Unhandled exception:\n{traceback.format_exc()}",
            expected=case.expected,
            duration_ms=duration,
        )
        _write_result(out_dir, result, case)
        return result

    duration = (time.time() - start) * 1000
    ok = case.matches(actual, case.expected)
    result = CaseResult(
        case_id=case.id,
        category=case.category,
        status=Status.PASS if ok else Status.FAIL,
        detail="matched expected outcome" if ok else "did NOT match expected outcome",
        actual=actual,
        expected=case.expected,
        duration_ms=duration,
    )
    _write_result(out_dir, result, case)
    return result


def run_all(cases: list) -> list:
    return [run_case(c) for c in cases]


def summarize(results: list) -> dict:
    counts = {s: 0 for s in Status}
    for r in results:
        counts[r.status] += 1
    return {s.value: n for s, n in counts.items()}
