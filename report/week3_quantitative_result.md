# Quantitative convergence result for the report

## Numerical verification paragraph

We addressed the review finding that the convergence evidence did not identify its theoretical benchmark, fitted uncertainty, or floating-point limit. For the constant-speed PMSM current subsystem, classical RK4 was run to the common final time (T=0.03\,\mathrm{s}) on uniform grids with (h=T/N). The reference value was evaluated independently from the affine system's matrix-exponential solution. We measured the componentwise endpoint error

\[
E(h)=\max_{j\in\{d,q\}} |i_{j,h}(T)-i_{j,\mathrm{ref}}(T)|.
\]

A least-squares fit of \(\log E\) against \(\log h\) over \(120\le N\le7680\) gave an observed order of \(4.038\pm0.009\), where the reported uncertainty is the regression standard error of the slope. This agrees with the theoretical fourth order. The finest points were excluded from the order fit because the error became non-monotone for \(h\lesssim9.77\times10^{-7}\,\mathrm{s}\), indicating a floating-point-dominated regime. The shaded 95% regression band is a diagnostic of the deterministic log-log fit; it is not repeated-sample statistical uncertainty.

## Figure caption

**Figure X. RK4 achieves fourth-order convergence before the floating-point floor.** Endpoint (L_\infty) current error is plotted against uniform step size for the illustrative constant-speed PMSM case. The fitted slope is (4.038\pm0.009) on the stated asymptotic range, consistent with the black order-4 reference line. The shaded fine-step region marks non-monotone floating-point behavior and is excluded from the fit. This matters because a fitted slope is credible only when its fitting range and numerical floor are visible.

![RK4 quantitative convergence figure](../figures/week3_rk4_quantitative.png)

## Scope sentence

This experiment verifies the implementation and observed order for the stated affine PMSM test case; it does not establish fourth-order behavior for nonsmooth switching inputs, state-dependent saturation models, or unresolved final motor parameters.
