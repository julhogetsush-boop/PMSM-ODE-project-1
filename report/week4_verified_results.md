# Week 4 reproduction and cost evidence

For the documented illustrative PMSM case, the RK4 endpoint-order estimate was recomputed as 4.03755492046426, matching the report to three significant digits (4.04). The Euler stability boundary is 0.2857142857142859 ms (0.286 ms to three significant digits). Reproduction uses the exact parameter file and the scripts shipped with this folder. The source SHA256 and environment are recorded with the run; a genuine Git commit must additionally be recorded before claiming the exact-commit Tutorial 4 requirement.

## Matched accuracy study

Cost is counted in right-hand-side evaluations, including implicit Newton residual evaluations. The same error ceiling is imposed on all methods; the first qualifying grid in the recorded powers-of-two sweep is selected. Errors are maximum component errors on each method's own time grid versus a Radau reference (rtol=1e-12, atol=1e-14). This is sampled trajectory error, not a proven supremum between nodes. Jacobian counts and Newton updates are reported separately. Jacobian assembly, linear solves, reference construction, setup and plotting are excluded from the primary count, so this experiment does not rank wall-clock time.

| Error ceiling [A] | Method | N | Measured error [A] | RHS calls |
|---:|---|---:|---:|---:|
| 0.01 | euler | 15360 | 0.00827688 | 15360 |
| 0.005 | euler | 30720 | 0.00413195 | 30720 |
| 0.01 | rk4 | 60 | 0.00473035 | 240 |
| 0.005 | rk4 | 60 | 0.00473035 | 240 |
| 0.01 | implicit | 15360 | 0.00822511 | 30720 |
| 0.005 | implicit | 30720 | 0.004119 | 61440 |

**Figure caption.** Cost and accuracy for the same constant-speed PMSM IVP. Each curve shows the measured error against Radau versus counted RHS evaluations. Dashed lines indicate the shared accuracy ceilings. RK4 reaches both ceilings with fewer RHS evaluations than either first-order method in this experiment, which matters because an efficiency comparison must control the attained accuracy.

![Cost and accuracy](../figures/week4_cost_accuracy.png)

The model has a pair of equally damped complex eigenvalues; its explicit-step restriction is primarily oscillatory, rather than evidence of widely separated decay scales. The cost result is specific to this affine teaching case.
