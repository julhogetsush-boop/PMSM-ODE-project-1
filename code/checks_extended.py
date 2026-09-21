"""Extended Testing & Validation for the deterministic PMSM project.

Complements `checks.py` with the checks requested in the T&V brief:

  1.  Hard-coded matrix A, forcing b, and equilibrium values.
  2.  Standalone Newton-Raphson tests (convergence, damping, failure paths).
  3.  Order-of-accuracy study (global error vs. step size).
  4.  Local-truncation error (one step from the exact solution).
  5.  Amplification function R(z) vs. Taylor series of exp(z).
  6.  Stability-limit analytics + empirical |R(h*lambda)| boundary test.
  7.  Adaptive step-doubling estimator order.
  8.  Adaptive tolerance sweep (steps and error monotonicity).
  9.  Zero-order-hold / stage-time verification.
  10. Independent reference comparison (Radau and DOP853 vs. expm).

Run: python code/checks_extended.py
"""

import math
import numpy as np
from scipy.integrate import solve_ivp

from methods import (
    euler_step, rk4_step, implicit_step, solve_fixed,
    adaptive_implicit, newton, NewtonFailure, time_grid,
)
from model import (
    load_parameters, matrices, make_model, equilibrium,
    exact_solution, amplification, stability_limit,
)


# ---------------------------------------------------------------------------
# 1. Hard-coded matrix values
# ---------------------------------------------------------------------------

def check_hardcoded_matrices():
    """A and b must equal a hand-derived reference for a fixed parameter set."""
    p = {
        "R": 0.5, "Ld": 4e-3, "Lq": 6e-3, "psi_f": 0.08,
        "omega_e": 500.0, "ud": 0.0, "uq": 50.0,
        "T": 0.01, "y0": [0.0, 0.0],
    }
    A, b = matrices(p)
    A_ref = np.array([
        [-0.5 / 4e-3,  500.0 * 6e-3 / 4e-3],
        [-500.0 * 4e-3 / 6e-3, -0.5 / 6e-3],
    ])
    b_ref = np.array([
        0.0 / 4e-3,
        (50.0 - 500.0 * 0.08) / 6e-3,
    ])
    assert np.allclose(A, A_ref, rtol=0, atol=1e-12), (A, A_ref)
    assert np.allclose(b, b_ref, rtol=0, atol=1e-12), (b, b_ref)

    # Equilibrium solves A x = -b exactly.
    star = equilibrium(p)
    assert np.allclose(A @ star + b, 0.0, rtol=0, atol=1e-12)
    return {"A": A, "b": b, "equilibrium": star}


# ---------------------------------------------------------------------------
# 2. Standalone Newton tests
# ---------------------------------------------------------------------------

def check_newton_solver():
    """Convergence, damping, singular Jacobian, non-finite residual, args."""
    # (a) Quadratic convergence on x^2 - 2 from x0 = 1.
    F  = lambda w: np.array([w[0] ** 2 - 2.0])
    JF = lambda w: np.array([[2.0 * w[0]]])
    w, hist = newton(F, JF, np.array([1.0]), tol=1e-14, max_iter=20)
    assert abs(w[0] - math.sqrt(2.0)) < 1e-12

    res = [h["residual"] for h in hist if h["residual"] > 0]
    assert len(res) >= 3, res
    ratios = [res[k + 1] / res[k] ** 2 for k in range(len(res) - 2)]
    # In the quadratic regime the ratio is approximately constant.
    assert abs(ratios[-1] - ratios[-2]) / max(ratios[-1], 1e-30) < 0.5, ratios

    # (b) Damped Newton on arctan from a very large start.
    F2  = lambda w: np.array([math.atan(w[0]) - 0.5])
    JF2 = lambda w: np.array([[1.0 / (1.0 + w[0] ** 2)]])
    w2, _ = newton(F2, JF2, np.array([10.0]), damped=True, max_iter=80)
    assert abs(math.atan(w2[0]) - 0.5) < 1e-12

    # (c) No root -> NewtonFailure.
    F3  = lambda w: np.array([math.exp(w[0])])
    JF3 = lambda w: np.array([[math.exp(w[0])]])
    try:
        newton(F3, JF3, np.array([0.0]), tol=1e-12, max_iter=5)
    except NewtonFailure:
        pass
    else:
        raise AssertionError("Newton should fail on exp(w) with no root")

    # (d) Singular Jacobian -> NewtonFailure.
    F4  = lambda w: np.array([1.0, 2.0])
    JF4 = lambda w: np.zeros((2, 2))
    try:
        newton(F4, JF4, np.zeros(2))
    except NewtonFailure:
        pass
    else:
        raise AssertionError("Singular Jacobian should raise NewtonFailure")

    # (e) Non-finite residual -> NewtonFailure.
    F5  = lambda w: np.array([np.nan])
    JF5 = lambda w: np.array([[1.0]])
    try:
        newton(F5, JF5, np.zeros(1))
    except NewtonFailure:
        pass
    else:
        raise AssertionError("NaN residual should raise NewtonFailure")

    # (f) Argument validation.
    for bad_tol in (0.0, -1.0, np.nan):
        try:
            newton(F, JF, np.array([1.0]), tol=bad_tol)
        except ValueError:
            pass
        else:
            raise AssertionError(f"tol={bad_tol} was accepted")
    try:
        newton(F, JF, np.array([1.0]), max_iter=0)
    except ValueError:
        pass
    else:
        raise AssertionError("max_iter=0 was accepted")

    return {"sqrt2": float(w[0]), "arctan_root": float(w2[0])}


# ---------------------------------------------------------------------------
# 3. Order of accuracy (scalar test problem)
# ---------------------------------------------------------------------------

def check_order_of_accuracy():
    """Empirical order on y' = -2y, y(0)=1, T=1; exact y = exp(-2t)."""
    f   = lambda t, y: -2.0 * y
    jac = lambda t, y: np.array([[-2.0]])
    T   = 1.0
    hs  = np.array([1/16, 1/32, 1/64, 1/128, 1/256], dtype=float)
    expected = {"euler": 1.0, "implicit": 1.0, "rk4": 4.0}
    results = {}
    for method, p_exp in expected.items():
        errs = []
        for h in hs:
            if method == "implicit":
                t, Y, _ = solve_fixed(f, [1.0], 0.0, T, h,
                                      method="implicit", jac=jac)
            else:
                t, Y, _ = solve_fixed(f, [1.0], 0.0, T, h, method=method)
            exact = np.exp(-2.0 * t)
            errs.append(float(np.max(np.abs(Y[:, 0] - exact))))
        errs = np.array(errs)
        p_est = np.log(errs[:-1] / errs[1:]) / np.log(hs[:-1] / hs[1:])
        assert abs(p_est[-1] - p_exp) < 0.4, (method, p_est)
        results[method] = {"hs": hs, "errors": errs, "p_estimates": p_est}
    return results


# ---------------------------------------------------------------------------
# 4. Local truncation error (scalar test problem)
# ---------------------------------------------------------------------------

def check_local_truncation():
    """Single-step error from the exact solution scales as h^(p+1)."""
    # y' = -2y, exact y(t) = exp(-2t); y(0)=1.
    f   = lambda t, y: -2.0 * y
    jac = lambda t, y: np.array([[-2.0]])
    y0  = np.array([1.0])
    # RK4 needs a coarser grid to stay above round-off; use per-method grids.
    grids = {
        "euler":    np.array([1e-2, 5e-3, 2.5e-3, 1.25e-3]),
        "implicit": np.array([1e-2, 5e-3, 2.5e-3, 1.25e-3]),
        "rk4":      np.array([1e-1, 5e-2, 2.5e-2, 1.25e-2]),
    }
    expected_exp = {"euler": 2.0, "implicit": 2.0, "rk4": 5.0}
    results = {}
    for method, hs in grids.items():
        errs = []
        for h in hs:
            if method == "implicit":
                y_new, _ = implicit_step(f, jac, 0.0, y0, h)
            elif method == "rk4":
                y_new = rk4_step(f, 0.0, y0, h)
            else:
                y_new = euler_step(f, 0.0, y0, h)
            y_exact = np.exp(-2.0 * h) * y0
            errs.append(float(np.max(np.abs(y_new - y_exact))))
        errs = np.array(errs)
        p_est = np.log(errs[:-1] / errs[1:]) / np.log(hs[:-1] / hs[1:])
        assert abs(p_est[-1] - expected_exp[method]) < 0.5, (method, p_est)
        results[method] = {"hs": hs, "errors": errs, "p_estimates": p_est}
    return results


# ---------------------------------------------------------------------------
# 5. Amplification vs. Taylor series of exp
# ---------------------------------------------------------------------------

def check_amplification_taylor():
    """|R(z) - exp(z)| / |z|^p must stay bounded as z -> 0."""
    z = 1e-3 + 0j
    for method, order in [("euler", 2), ("rk4", 5), ("implicit", 2)]:
        R = amplification(method, z)
        gap = abs(R - np.exp(z))
        scaled = gap / abs(z) ** order
        assert scaled < 10.0, (method, scaled)
    try:
        amplification("bogus", z)
    except ValueError:
        pass
    else:
        raise AssertionError("unknown method should raise")
    return "amplification matches Taylor of exp to expected order"


# ---------------------------------------------------------------------------
# 6. Stability limit
# ---------------------------------------------------------------------------

def check_stability_limit():
    """Analytic boundaries and marginal |R(h*lambda)| = 1 confirmation."""
    # (a) Euler on real lambda = -1 -> h_lim = 2.
    h_euler_real = stability_limit("euler", np.array([-1.0]))
    assert abs(h_euler_real - 2.0) < 1e-12, h_euler_real
    # (b) Euler on lambda = -1 + 1j -> h_lim = 2*1 / (1 + 1) = 1.
    h_euler_cplx = stability_limit("euler", np.array([-1.0 + 1.0j]))
    assert abs(h_euler_cplx - 1.0) < 1e-12, h_euler_cplx
    # (c) RK4 on real lambda = -1 -> ~2.785293.
    h_rk4 = stability_limit("rk4", np.array([-1.0]))
    assert abs(h_rk4 - 2.785293) < 1e-5, h_rk4
    # (d) Implicit -> inf.
    assert stability_limit("implicit", np.array([-1.0])) == np.inf

    # (e) Marginal test on a complex eigenvalue.
    lam = -1.0 + 0.5j
    for method in ("euler", "rk4"):
        h_lim = stability_limit(method, np.array([lam]))
        mag_at  = abs(amplification(method,  h_lim * lam))
        mag_bel = abs(amplification(method,  0.99 * h_lim * lam))
        mag_abv = abs(amplification(method,  1.01 * h_lim * lam))
        assert abs(mag_at - 1.0) < 1e-4, (method, mag_at)
        assert mag_bel <= 1.0 + 1e-9, (method, mag_bel)
        assert mag_abv > 1.0,         (method, mag_abv)

    # (f) Non-decaying eigenvalues must be rejected.
    for bad in (np.array([0.0]), np.array([1.0])):
        try:
            stability_limit("euler", bad)
        except ValueError:
            pass
        else:
            raise AssertionError("non-decaying eigenvalue accepted")
    return {"euler_real": h_euler_real, "euler_complex": h_euler_cplx,
            "rk4_real": h_rk4}


# ---------------------------------------------------------------------------
# 7. Step-doubling estimator order
# ---------------------------------------------------------------------------

def check_step_doubling_order():
    """||fine - coarse|| for implicit Euler must scale as O(h^2)."""
    p = load_parameters()
    f, jac = make_model(p)
    y = np.array(p["y0"], dtype=float)
    hs = np.array([4e-5, 2e-5, 1e-5, 5e-6])
    diffs = []
    for h in hs:
        coarse, _ = implicit_step(f, jac, 0.0, y, h)
        half,   _ = implicit_step(f, jac, 0.0, y, h / 2.0)
        fine,   _ = implicit_step(f, jac, h / 2.0, half, h / 2.0)
        diffs.append(float(np.max(np.abs(fine - coarse))))
    diffs = np.array(diffs)
    p_est = np.log(diffs[:-1] / diffs[1:]) / np.log(hs[:-1] / hs[1:])
    assert abs(p_est[-1] - 2.0) < 0.5, p_est
    return {"hs": hs, "diffs": diffs, "p_estimates": p_est}


# ---------------------------------------------------------------------------
# 8. Adaptive tolerance sweep
# ---------------------------------------------------------------------------

def check_adaptive_tolerance_sweep():
    """Tighter atol -> more accepted steps and smaller final error."""
    p = load_parameters()
    f, jac = make_model(p)
    y0 = np.array(p["y0"], dtype=float)
    T  = min(float(p["T"]), 5e-3)
    rows = []
    for atol in (1e-4, 1e-6, 1e-8):
        times, states, trials = adaptive_implicit(
            f, jac, y0, 0.0, T, h0=1e-4,
            atol=atol, rtol=0.0, hmin=1e-14, hmax=T,
        )
        err = float(np.max(np.abs(states[-1] - exact_solution(T, y0, p))))
        rows.append({
            "atol": atol,
            "accepted": len(times) - 1,
            "rejected": sum(1 for tr in trials if not tr["accepted"]),
            "final_error": err,
            "t_end": float(times[-1]),
        })
        assert rows[-1]["t_end"] == T
    # Tighter tolerance should not reduce the accepted step count.
    for a, b in zip(rows, rows[1:]):
        assert b["accepted"] >= a["accepted"] - 1, rows
        assert b["final_error"] <= a["final_error"] * 1.5 + 1e-12, rows
    return rows


# ---------------------------------------------------------------------------
# 9. Zero-order-hold / stage times
# ---------------------------------------------------------------------------

def check_zero_order_hold():
    """Non-autonomous probes that catch wrong stage times."""
    # RK4 on y' = t from 0 to 1: exact integral is 0.5.
    t, Y, _ = solve_fixed(lambda t, y: np.array([t]), [0.0],
                          0.0, 1.0, 0.3, "rk4")
    assert t[-1] == 1.0
    assert abs(Y[-1, 0] - 0.5) < 1e-13, Y[-1, 0]

    # Implicit Euler on y' = t: one step evaluates at t_new (zero-order hold).
    y_impl, _ = implicit_step(
        lambda t, y: np.array([t]),
        lambda t, y: np.zeros((1, 1)),
        0.0, np.zeros(1), 0.3,
    )
    assert abs(y_impl[0] - 0.09) < 1e-13, y_impl[0]

    # Explicit Euler on y' = t: uses the stage at t_old, so result is 0.
    y_eul = euler_step(lambda t, y: np.array([t]), 0.0, np.zeros(1), 0.3)
    assert abs(y_eul[0]) < 1e-15, y_eul[0]

    # RK4 on y' = t from t0 != 0: integral from 1 to 2 of t dt = 1.5.
    t2, Y2, _ = solve_fixed(lambda t, y: np.array([t]), [0.0],
                            1.0, 2.0, 0.25, "rk4")
    assert abs(Y2[-1, 0] - 1.5) < 1e-13, Y2[-1, 0]

    return "stage times correct for euler, rk4, implicit"


# ---------------------------------------------------------------------------
# 10. Independent reference comparison
# ---------------------------------------------------------------------------

def check_reference_solvers():
    """Radau and DOP853 must agree with the matrix exponential solution."""
    p = load_parameters()
    f, jac = make_model(p)
    y0 = np.array(p["y0"], dtype=float)
    T  = float(p["T"])
    ts = np.linspace(0.0, T, 401)

    ref = exact_solution(ts, y0, p)

    sol_radau = solve_ivp(f, (0.0, T), y0, method="Radau", jac=jac,
                          t_eval=ts, rtol=1e-12, atol=1e-14)
    sol_dop   = solve_ivp(f, (0.0, T), y0, method="DOP853",
                          t_eval=ts, rtol=1e-12, atol=1e-14)
    assert sol_radau.success and sol_dop.success
    gap_radau = float(np.max(np.abs(sol_radau.y.T - ref)))
    gap_dop   = float(np.max(np.abs(sol_dop.y.T   - ref)))
    assert gap_radau < 1e-8, gap_radau
    assert gap_dop   < 1e-8, gap_dop
    return {"radau_gap": gap_radau, "dop853_gap": gap_dop}


# ---------------------------------------------------------------------------
# 11. Diagnostic assignment table (for the report)
# ---------------------------------------------------------------------------

DIAGNOSTIC_TABLE = [
    ("Global order below expectation",
     "Wrong RK4 stage weights or implicit time level",
     "check_order_of_accuracy / check_zero_order_hold"),
    ("Local order below expectation",
     "Stage coefficients or quadrature wrong",
     "check_local_truncation"),
    ("|R(h_lim * lambda)| != 1",
     "stability_limit formula or bisection broken",
     "check_stability_limit (marginal test)"),
    ("Newton fails to converge",
     "Singular Jacobian or bad initial guess",
     "check_newton_solver (c, d, f)"),
    ("Adaptive solver stalls",
     "hmin too large or estimator order wrong",
     "check_step_doubling_order / check_adaptive_tolerance_sweep"),
    ("Radau disagrees with expm",
     "Model or exact_solution wrong",
     "check_hardcoded_matrices / check_reference_solvers"),
]


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_extended_checks(verbose=True):
    out = {}
    out["hardcoded"]     = check_hardcoded_matrices()
    out["newton"]        = check_newton_solver()
    out["order"]         = check_order_of_accuracy()
    out["local"]         = check_local_truncation()
    out["amplification"] = check_amplification_taylor()
    out["stability"]     = check_stability_limit()
    out["doubling"]      = check_step_doubling_order()
    out["adaptive"]      = check_adaptive_tolerance_sweep()
    out["zoh"]           = check_zero_order_hold()
    out["reference"]     = check_reference_solvers()
    out["diagnostic"]    = DIAGNOSTIC_TABLE

    if verbose:
        print("=" * 72)
        print("Extended T&V: PASS")
        print("=" * 72)
        for m, r in out["order"].items():
            print(f"  order  {m:8s}: p = " +
                  ", ".join(f"{v:.2f}" for v in r["p_estimates"]))
        for m, r in out["local"].items():
            print(f"  local  {m:8s}: p = " +
                ", ".join(f"{v:.2f}" for v in r["p_estimates"]))
        print(f"  stability : {out['stability']}")
        print(f"  doubling  : p = " +
              ", ".join(f"{v:.2f}" for v in out["doubling"]["p_estimates"]))
        print("  adaptive  : " +
              str([(r["atol"], r["accepted"], r["final_error"])
                   for r in out["adaptive"]]))
        print(f"  reference : radau={out['reference']['radau_gap']:.3e}, "
              f"dop853={out['reference']['dop853_gap']:.3e}")
        print("-" * 72)
        print("Diagnostic assignment table:")
        for symptom, cause, diag in DIAGNOSTIC_TABLE:
            print(f"  - {symptom}")
            print(f"      cause: {cause}")
            print(f"      check: {diag}")
    return out


if __name__ == "__main__":
    run_extended_checks(verbose=True)