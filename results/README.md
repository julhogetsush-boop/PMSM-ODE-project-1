# Reproducible results

Run `python code/run_all.py` from the repository root. This directory then contains:

- `summary.json`: all summary values, illustrative PMSM parameters, validation status and actual Python/package versions.
- `pmsm_convergence.csv`: all five grids for each of the three methods, whole-grid and endpoint errors.
- `pmsm_stability_sweep.csv`: the actual step sizes, amplification radii and perturbation ratios.
- `pmsm_adaptive_trials.csv`: every accepted/rejected trial for the middle tolerance.
- `pmsm_adaptive_summary.csv`: three tolerances, global errors and RHS counts.
- `pmsm_newton_trace.csv`: the constant-coefficient model's first-step Newton residual.
- `synthetic_stiff_case.csv`: the separately labelled stiff parameter experiment.
- `lab_*.csv`: the classroom exercises, RC first-step histories and full step-size scan.

Raw generated CSV/JSON files are included in the delivered local package but ignored by Git. The generated, readable result note and selected final PNG figures are tracked for browsing. Every value can be regenerated; none is a hard-coded copied benchmark.
