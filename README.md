# PMSM ODE Project 1 - Algorithm Implementation

This repository folder contains the report-aligned numerical implementation for
the PMSM electrical initial value problem.

The implemented methods are:

- Explicit Euler
- Classical RK4
- Implicit Euler solved by Newton iteration
- A linear implicit Euler cross-check for the affine PMSM benchmark
- Adaptive RK4 with step doubling and componentwise RMS error scaling

The benchmark model is

```text
i'(t) = A i(t) + b,        i = [i_d, i_q]^T,
```

using the parameters declared in the report:

```text
Rs = 0.5 ohm, Ld = 0.004 H, Lq = 0.006 H,
psi_f = 0.08 Wb, omega_e = 500 rad/s,
(u_d, u_q) = (0, 50) V, i(0) = (0, 0) A, T = 0.05 s.
```

## How to Run

From the repository root:

```powershell
python code/run_all.py
```

On Hannah's computer, if the default `python` points to the unstable alpha
interpreter, use:

```powershell
C:\Users\Yolanda ZHANG\Documents\Codex\2026-09-08\new-chat\work\tools\Python313\python.exe code\run_all.py
```

The script runs 13 numerical audit groups and regenerates:

- `results/summary.json`
- `results/checks.json`
- `experiment_manifest.json`
- 24 CSV tables in `results/`
- 14 PNG figures in `figures/`
- reproduction notes in `review/`

## Individual Challenge Candidate

A reproducible line from the current code is:

```text
Results from your run: Explicit Euler observed order = 1.0163802167...; standard error = 0.00344105046.... This rounds to 1.016380.
```

Use the final GitHub commit ID after uploading this folder to the team
repository.

## Scope

This package supports the Algorithm implementation role. It does not claim new
hardware validation or a rerun of the historical Simulink model. Synthetic
interface-mismatch experiments are labelled as synthetic in the generated
results and review notes.
