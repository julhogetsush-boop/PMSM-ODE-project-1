"""One-command regeneration entry point for the report-aligned package."""

from __future__ import annotations

import json
import sys
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import checks
import experiments


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    (PROJECT_ROOT / "results").mkdir(exist_ok=True)

    check_records = checks.write_checks(PROJECT_ROOT)
    failed = [record for record in check_records if not record["passed"]]
    if failed:
        for record in failed:
            print(f"[FAIL] {record['name']}: {record.get('value')}")
        raise SystemExit("Numerical audit checks failed.")

    print(f"Passed {len(check_records)} numerical audit groups.")
    summary = experiments.run_all(PROJECT_ROOT)
    summary["checks"] = check_records
    (PROJECT_ROOT / "results" / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    explicit = summary["convergence"]["fits"]["explicit_euler"]
    rk4_report = summary["report_table_3_values"]["rk4"]
    implicit = summary["convergence"]["fits"]["implicit_euler"]
    print("PMSM Algorithm implementation experiments complete.")
    print(f"Explicit Euler observed order = {explicit['order']}; standard error = {explicit['slope_standard_error']}.")
    print(f"RK4 observed order = {rk4_report['order']}; standard error = {rk4_report['slope_standard_error']}.")
    print(f"Implicit Euler observed order = {implicit['order']}; standard error = {implicit['slope_standard_error']}.")
    print(f"Explicit Euler h_max = {summary['stability']['h_explicit']:.6e}.")
    print(f"RK4 h_max = {summary['stability']['h_rk4']:.6e}.")
    print("Wrote results/summary.json, experiment_manifest.json, 24 CSV files and 14 figures.")


if __name__ == "__main__":
    main()
