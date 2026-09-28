"""Numerical audit checks for the PMSM implementation."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np

from experiments import (
    T_SPAN,
    Y0,
    explicit_euler_hmax,
    max_mesh_error,
    rk4_hmax,
)
from model import (
    DEFAULT_PARAMETERS,
    PMSMParameters,
    equilibrium,
    exact_solution,
    flux_error_norm,
    forcing_vector,
    jacobian,
    rhs,
)
from solvers import (
    adaptive_rk4,
    explicit_euler,
    implicit_euler_linear,
    implicit_euler_newton,
    rk4,
)


def _record(name: str, passed: bool, value: Any = None, limit: Any = None) -> Dict[str, Any]:
    return {"name": name, "passed": bool(passed), "value": value, "limit": limit}


def run_checks(project_root: Path | None = None) -> List[Dict[str, Any]]:
    """Run 13 audit groups and return machine-readable check records."""
    p = DEFAULT_PARAMETERS
    A = jacobian(p)
    b = forcing_vector(0.0, p)
    f = lambda t, y: rhs(t, y, p)
    jac = lambda t, y: A
    checks: List[Dict[str, Any]] = []

    # 1. Matrix and forcing match the declared benchmark.
    checks.append(_record("model_matrix_shape", A.shape == (2, 2) and b.shape == (2,), {"A": A.tolist(), "b": b.tolist()}))

    # 2. Exact solution respects the initial condition and equilibrium.
    init_err = float(np.linalg.norm(exact_solution(0.0, Y0, p) - Y0))
    eq_res = float(np.linalg.norm(A @ equilibrium(p) + b))
    checks.append(_record("exact_initial_and_equilibrium", init_err < 1e-14 and eq_res < 1e-12, {"initial_error": init_err, "equilibrium_residual": eq_res}, 1e-12))

    # 3. Zero-speed case reduces to two independent RL equations.
    p0 = PMSMParameters(omega_e=0.0)
    t = np.linspace(0.0, 0.05, 101)
    y_exact = exact_solution(t, Y0, p0)
    rl = np.column_stack(
        [
            (p0.u_d / p0.Rs) * (1 - np.exp(-(p0.Rs / p0.Ld) * t)),
            (p0.u_q / p0.Rs) * (1 - np.exp(-(p0.Rs / p0.Lq) * t)),
        ]
    )
    rl_err = float(np.max(np.linalg.norm(y_exact - rl, axis=1)))
    checks.append(_record("omega_zero_RL_reference", rl_err < 1e-12, rl_err, 1e-12))

    # 4. Explicit Euler is first order on the declared fitted range.
    e_h = np.array([50e-6, 25e-6, 12.5e-6, 6.25e-6])
    e_err = np.array([max_mesh_error(explicit_euler(f, Y0, T_SPAN, float(h)), p) for h in e_h])
    e_order = float(np.polyfit(np.log(e_h), np.log(e_err), 1)[0])
    checks.append(_record("explicit_euler_order", abs(e_order - 1.0) < 0.05, e_order, "within 0.05 of 1"))

    # 5. RK4 is fourth order before roundoff dominates.
    r_err = np.array([max_mesh_error(rk4(f, Y0, T_SPAN, float(h)), p) for h in e_h])
    r_order = float(np.polyfit(np.log(e_h), np.log(r_err), 1)[0])
    checks.append(_record("rk4_order", abs(r_order - 4.0) < 0.05, r_order, "within 0.05 of 4"))

    # 6. Implicit Euler is first order.
    i_err = np.array([max_mesh_error(implicit_euler_linear(A, b, Y0, T_SPAN, float(h)), p) for h in e_h])
    i_order = float(np.polyfit(np.log(e_h), np.log(i_err), 1)[0])
    checks.append(_record("implicit_euler_order", abs(i_order - 1.0) < 0.05, i_order, "within 0.05 of 1"))

    # 7. Newton implicit Euler matches the affine linear solve.
    newton = implicit_euler_newton(f, jac, Y0, (0.0, 0.005), 1e-4)
    linear = implicit_euler_linear(A, b, Y0, (0.0, 0.005), 1e-4)
    diff = float(np.max(np.linalg.norm(newton.y - linear.y, axis=1)))
    checks.append(_record("implicit_newton_matches_linear_affine_solve", diff < 1e-11, diff, 1e-11))

    # 8. Euler and RK4 stability limits are positive and ordered as expected.
    h_e = explicit_euler_hmax(A)
    h_r = rk4_hmax(A)
    checks.append(_record("stability_step_limits", 0.0007 < h_e < 0.0009 and h_r > h_e, {"euler": h_e, "rk4": h_r}))

    # 9. Deliberately too-large Euler step grows relative to a stable step.
    stable = explicit_euler(f, Y0, T_SPAN, 0.5 * h_e)
    unstable = explicit_euler(f, Y0, T_SPAN, 1.25 * h_e)
    checks.append(_record("explicit_euler_deliberate_instability", np.linalg.norm(unstable.y[-1]) > np.linalg.norm(stable.y[-1]), {"stable_final_norm": float(np.linalg.norm(stable.y[-1])), "unstable_final_norm": float(np.linalg.norm(unstable.y[-1]))}))

    # 10. Flux error norm contracts for the exact trajectory toward equilibrium.
    t_flux = np.linspace(0.0, 0.05, 200)
    flux = flux_error_norm(exact_solution(t_flux, Y0, p), equilibrium(p), p)
    checks.append(_record("flux_norm_contracts", bool(np.all(np.diff(flux) <= 1e-14)), float(np.max(np.diff(flux))), "<= 0"))

    # 11. Adaptive RK4 completes with accepted and rejected trials under target.
    adaptive = adaptive_rk4(f, Y0, T_SPAN, 0.005, atol=1e-5, rtol=1e-5, h_max=0.005)
    adaptive_err = max_mesh_error(adaptive, p)
    checks.append(_record("adaptive_rk4_global_check", adaptive_err < 1e-3 and adaptive.stats.steps > 0, {"Emax": adaptive_err, "accepted": adaptive.stats.steps, "rejected": adaptive.stats.rejected_steps}, 1e-3))

    # 12. Input-bias mismatch persists even as RK4 numerical error shrinks.
    p_unbiased = PMSMParameters(Rs=0.36, Ld=0.0002, Lq=0.0002, psi_f=0.00639541518689, omega_e=1500.0, u_d=0.0, u_q=12.0)
    p_biased = PMSMParameters(Rs=0.36, Ld=0.0002, Lq=0.0002, psi_f=0.00639541518689, omega_e=1500.0, u_d=0.1, u_q=12.0)
    res = rk4(lambda t, y: rhs(t, y, p_biased), Y0, (0.0, 0.02), 10e-6)
    biased_ref = exact_solution(res.t, Y0, p_biased)
    unbiased_ref = exact_solution(res.t, Y0, p_unbiased)
    numerical = float(np.max(np.linalg.norm(res.y - biased_ref, axis=1)))
    mismatch = float(np.linalg.norm(res.y[-1] - unbiased_ref[-1]))
    checks.append(_record("synthetic_input_bias_separates_error_sources", numerical < 1e-7 and mismatch > 0.1, {"numerical": numerical, "mismatch": mismatch}))

    # 13. Unequal-scale toy problem demonstrates scaling policy risk.
    def toy_rhs(t: float, y: np.ndarray) -> np.ndarray:
        return np.array([0.0, 100.0 * np.cos(100.0 * t)])

    componentwise = adaptive_rk4(toy_rhs, np.array([1e6, 0.0]), (0.0, 0.2), 0.1, atol=np.array([1e-6, 1e-6]), rtol=1e-6, h_max=0.1)
    shared = adaptive_rk4(toy_rhs, np.array([1e6, 0.0]), (0.0, 0.2), 0.1, atol=1.000001, rtol=1e-15, h_max=0.1)
    c_err = float(np.max(np.abs(componentwise.y[:, 1] - np.sin(100.0 * componentwise.t))))
    s_err = float(np.max(np.abs(shared.y[:, 1] - np.sin(100.0 * shared.t))))
    checks.append(_record("unequal_scale_counterexample", s_err > 0.1 and c_err < 1e-4, {"componentwise": c_err, "shared": s_err}))

    return checks


def write_checks(project_root: Path) -> List[Dict[str, Any]]:
    checks = run_checks(project_root)
    out = project_root / "results" / "checks.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(__import__("json").dumps(checks, indent=2), encoding="utf-8")
    return checks
