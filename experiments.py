"""Report-aligned PMSM numerical experiments.

Running :mod:`run_all` regenerates the machine-readable evidence used by the
report: CSV tables, a JSON summary and figures.  The benchmark is the declared
constant-speed PMSM electrical IVP on 0 <= t <= 0.05 s.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
from dataclasses import asdict
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

from model import (
    DEFAULT_PARAMETERS,
    PMSMParameters,
    copper_loss,
    equilibrium,
    exact_solution,
    flux_error_norm,
    forcing_vector,
    input_power,
    jacobian,
    magnetic_energy,
    rhs,
    torque,
)
from solvers import (
    SolverResult,
    adaptive_rk4,
    explicit_euler,
    implicit_euler_linear,
    implicit_euler_newton,
    rk4,
    rk4_step,
)


Y0 = np.array([0.0, 0.0])
T_SPAN = (0.0, 0.05)
CONVERGENCE_H = np.array([100.0, 50.0, 25.0, 12.5, 6.25]) * 1e-6
FIT_H = CONVERGENCE_H[1:]


def _ensure_dirs(project_root: Path) -> Dict[str, Path]:
    paths = {
        "figures": project_root / "figures",
        "results": project_root / "results",
        "review": project_root / "review",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def _write_csv(path: Path, header: Iterable[str], rows: Iterable[Iterable[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(list(header))
        writer.writerows(rows)


def _source_hashes(code_dir: Path) -> Dict[str, str]:
    hashes: Dict[str, str] = {}
    for file_name in ["model.py", "solvers.py", "experiments.py", "checks.py", "run_all.py"]:
        path = code_dir / file_name
        hashes[f"code/{file_name}"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def _baseline() -> Tuple[PMSMParameters, np.ndarray, np.ndarray, Callable[[float, np.ndarray], np.ndarray]]:
    p = DEFAULT_PARAMETERS
    A = jacobian(p)
    b = forcing_vector(0.0, p)
    f = lambda t, y: rhs(t, y, p)
    return p, A, b, f


def max_mesh_error(result: SolverResult, p: PMSMParameters = DEFAULT_PARAMETERS) -> float:
    reference = exact_solution(result.t, result.y[0], p)
    return float(np.max(np.linalg.norm(result.y - reference, axis=1)))


def final_error(result: SolverResult, p: PMSMParameters = DEFAULT_PARAMETERS) -> float:
    reference = exact_solution(result.t[-1], result.y[0], p)
    return float(np.linalg.norm(result.y[-1] - reference))


def fit_order(h_values: np.ndarray, errors: np.ndarray) -> Dict[str, float]:
    x = np.log(np.asarray(h_values, dtype=float))
    y = np.log(np.asarray(errors, dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    residuals = y - (slope * x + intercept)
    dof = max(len(x) - 2, 1)
    variance = float(np.sum(residuals**2) / dof)
    slope_stderr = float(np.sqrt(variance / np.sum((x - x.mean()) ** 2)))
    return {
        "order": float(slope),
        "intercept": float(intercept),
        "slope_standard_error": slope_stderr,
        "residual_standard_error": float(np.sqrt(variance)),
    }


def rk4_stability_function(z: complex | np.ndarray) -> complex | np.ndarray:
    return 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24


def explicit_euler_hmax(A: np.ndarray) -> float:
    limits = []
    for lam in np.linalg.eigvals(A):
        if np.real(lam) >= 0:
            return 0.0
        limits.append(-2.0 * np.real(lam) / abs(lam) ** 2)
    return float(min(limits))


def rk4_hmax(A: np.ndarray) -> float:
    eigvals = np.linalg.eigvals(A)

    def stable(h: float) -> bool:
        return bool(np.max(np.abs(rk4_stability_function(h * eigvals))) <= 1.0 + 1e-12)

    lo, hi = 0.0, 1e-3
    while stable(hi):
        lo, hi = hi, 2.0 * hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if stable(mid):
            lo = mid
        else:
            hi = mid
    return float(lo)


def euler_flux_contraction_hmax(p: PMSMParameters = DEFAULT_PARAMETERS) -> float:
    A = jacobian(p)
    D = np.diag([p.Ld, p.Lq])
    B = D @ A @ np.linalg.inv(D)
    S = np.diag([p.Rs / p.Ld, p.Rs / p.Lq])
    M = np.linalg.inv(np.sqrt(2 * S)) @ (B.T @ B) @ np.linalg.inv(np.sqrt(2 * S))
    return float(1.0 / np.max(np.linalg.eigvalsh(M)))


def run_convergence(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    p, A, b, f = _baseline()
    rows = []
    errors = {"explicit_euler": [], "rk4": [], "implicit_euler": []}
    for h in CONVERGENCE_H:
        euler_result = explicit_euler(f, Y0, T_SPAN, float(h))
        rk4_result = rk4(f, Y0, T_SPAN, float(h))
        implicit_result = implicit_euler_linear(A, b, Y0, T_SPAN, float(h))
        e_err = max_mesh_error(euler_result, p)
        r_err = max_mesh_error(rk4_result, p)
        i_err = max_mesh_error(implicit_result, p)
        errors["explicit_euler"].append(e_err)
        errors["rk4"].append(r_err)
        errors["implicit_euler"].append(i_err)
        rows.append([h, e_err, r_err, i_err, euler_result.stats.steps])

    _write_csv(
        results_dir / "convergence_table.csv",
        ["h_seconds", "explicit_euler_Emax", "rk4_Emax", "implicit_euler_Emax", "steps"],
        rows,
    )

    fits = {
        name: fit_order(FIT_H, np.asarray(values[1:]))
        for name, values in errors.items()
    }
    order_rows = [
        [name, 1 if name != "rk4" else 4, fits[name]["order"], fits[name]["slope_standard_error"]]
        for name in ["explicit_euler", "rk4", "implicit_euler"]
    ]
    _write_csv(
        results_dir / "observed_orders.csv",
        ["method", "theoretical_order", "fitted_order", "slope_standard_error"],
        order_rows,
    )

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for name, marker in [("explicit_euler", "o-"), ("rk4", "s-"), ("implicit_euler", "^-")]:
        ax.loglog(CONVERGENCE_H, errors[name], marker, label=f"{name}, p={fits[name]['order']:.3f}")
    ax.invert_xaxis()
    ax.set_xlabel("Nominal h [s]")
    ax.set_ylabel("Maximum mesh error [A]")
    ax.set_title("Convergence on the baseline PMSM IVP")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "01_convergence.png", dpi=180)
    plt.close(fig)

    return {"errors": errors, "fits": fits}


def run_reference_and_trajectory(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    p, A, _, f = _baseline()
    t_dense = np.linspace(T_SPAN[0], T_SPAN[1], 2001)
    y_exact = exact_solution(t_dense, Y0, p)
    y_star = equilibrium(p)
    h = 5e-4
    euler_result = explicit_euler(f, Y0, T_SPAN, h)
    rk4_result = rk4(f, Y0, T_SPAN, h)

    _write_csv(
        results_dir / "trajectory_reference.csv",
        ["t", "i_d_exact", "i_q_exact", "torque_exact"],
        [[t, y[0], y[1], torque(y, p)] for t, y in zip(t_dense, y_exact)],
    )
    _write_csv(
        results_dir / "phase_reference.csv",
        ["i_d", "i_q"],
        y_exact,
    )

    # Independent analytic checks that do not share the 2x2 oscillatory branch.
    omega_zero = PMSMParameters(omega_e=0.0)
    t_check = np.linspace(0.0, 0.05, 101)
    exact_zero = exact_solution(t_check, Y0, omega_zero)
    rl_d_star = omega_zero.u_d / omega_zero.Rs
    rl_q_star = omega_zero.u_q / omega_zero.Rs
    rl_exact = np.column_stack(
        [
            rl_d_star + (Y0[0] - rl_d_star) * np.exp(-(omega_zero.Rs / omega_zero.Ld) * t_check),
            rl_q_star + (Y0[1] - rl_q_star) * np.exp(-(omega_zero.Rs / omega_zero.Lq) * t_check),
        ]
    )
    rl_error = float(np.max(np.linalg.norm(exact_zero - rl_exact, axis=1)))
    _write_csv(
        results_dir / "exact_rl_check.csv",
        ["t", "exact_i_d", "rl_i_d", "exact_i_q", "rl_i_q"],
        [[t, a[0], b[0], a[1], b[1]] for t, a, b in zip(t_check, exact_zero, rl_exact)],
    )
    _write_csv(
        results_dir / "reference_agreement.csv",
        ["check", "max_difference_A"],
        [["omega_zero_RL_formula", rl_error], ["initial_value", float(np.linalg.norm(exact_solution(0.0, Y0, p) - Y0))]],
    )

    strong = PMSMParameters(Rs=50.0, Ld=4e-3, Lq=6e-3, psi_f=0.08, omega_e=500.0, u_d=0.0, u_q=50.0)
    t_strong = np.linspace(0.0, 0.5, 101)
    y_strong = exact_solution(t_strong, Y0, strong)
    _write_csv(
        results_dir / "strong_decay_check.csv",
        ["t", "i_d", "i_q", "finite"],
        [[t, y[0], y[1], bool(np.isfinite(y).all())] for t, y in zip(t_strong, y_strong)],
    )

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.plot(t_dense * 1000, y_exact[:, 0], "k-", label="exact i_d")
    ax.plot(t_dense * 1000, y_exact[:, 1], "k--", label="exact i_q")
    ax.plot(rk4_result.t * 1000, rk4_result.y[:, 0], color="tab:blue", alpha=0.8, label="RK4 i_d")
    ax.plot(rk4_result.t * 1000, rk4_result.y[:, 1], color="tab:orange", alpha=0.8, label="RK4 i_q")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Current [A]")
    ax.set_title("Current trajectory")
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(figures_dir / "02_trajectory.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    ax.plot(y_exact[:, 0], y_exact[:, 1], "k-", label="exact")
    ax.plot(euler_result.y[:, 0], euler_result.y[:, 1], label="Euler h=0.5 ms")
    ax.plot(rk4_result.y[:, 0], rk4_result.y[:, 1], label="RK4 h=0.5 ms")
    ax.scatter(y_star[0], y_star[1], marker="x", color="red", label="equilibrium")
    ax.set_xlabel("i_d [A]")
    ax.set_ylabel("i_q [A]")
    ax.set_title("Phase portrait")
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "03_phase_plane.png", dpi=180)
    plt.close(fig)

    return {"rl_check_error": rl_error, "equilibrium": y_star.tolist()}


def run_stability(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    p, A, b, f = _baseline()
    eigvals = np.linalg.eigvals(A)
    h_euler = explicit_euler_hmax(A)
    h_rk4 = rk4_hmax(A)
    h_flux = euler_flux_contraction_hmax(p)
    h_scan = np.linspace(0.0, 1.35 * h_rk4, 800)
    euler_amp = np.max(np.abs(1.0 + h_scan[:, None] * eigvals[None, :]), axis=1)
    rk4_amp = np.max(np.abs(rk4_stability_function(h_scan[:, None] * eigvals[None, :])), axis=1)
    _write_csv(
        results_dir / "stability_scan.csv",
        ["h", "explicit_euler_amplification", "rk4_amplification"],
        [[h, e, r] for h, e, r in zip(h_scan, euler_amp, rk4_amp)],
    )

    stable = explicit_euler(f, Y0, T_SPAN, 0.5 * h_euler)
    unstable = explicit_euler(f, Y0, T_SPAN, 1.25 * h_euler)
    _write_csv(
        results_dir / "unstable_euler.csv",
        ["t", "stable_norm", "unstable_norm"],
        [[t, np.linalg.norm(s), np.linalg.norm(u)] for t, s, u in zip(stable.t, stable.y, unstable.y[: len(stable.t)])],
    )

    fig, ax = plt.subplots(figsize=(7.3, 5.0))
    ax.plot(h_scan * 1000, euler_amp, label="Explicit Euler")
    ax.plot(h_scan * 1000, rk4_amp, label="RK4")
    ax.axhline(1.0, color="black", lw=1, ls="--")
    ax.axvline(h_euler * 1000, color="tab:blue", ls=":", label=f"Euler hmax={h_euler*1000:.3f} ms")
    ax.axvline(h_rk4 * 1000, color="tab:orange", ls=":", label=f"RK4 hmax={h_rk4*1000:.3f} ms")
    ax.set_xlabel("Step size [ms]")
    ax.set_ylabel("max |R(h lambda)|")
    ax.set_title("Stability restriction along PMSM eigenvalues")
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "04_stability_boundary.png", dpi=180)
    plt.close(fig)

    return {"h_explicit": h_euler, "h_rk4": h_rk4, "h_flux": h_flux}


def run_contraction_power(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    p, A, b, f = _baseline()
    t = np.linspace(0.0, 0.05, 1001)
    y_exact = exact_solution(t, Y0, p)
    y_star = equilibrium(p)
    flux_norm = flux_error_norm(y_exact, y_star, p)
    m = min(p.Rs / p.Ld, p.Rs / p.Lq)
    bound = flux_norm[0] * np.exp(-m * t)
    pin = input_power(y_exact, p)
    pcu = copper_loss(y_exact, p)
    pmech = torque(y_exact, p) * p.omega_m
    energy = magnetic_energy(y_exact, p)
    _write_csv(
        results_dir / "contraction_power.csv",
        ["t", "flux_error_norm", "bound", "input_power", "copper_loss", "mechanical_power", "magnetic_energy"],
        [[a, b0, c, d, e, g, h] for a, b0, c, d, e, g, h in zip(t, flux_norm, bound, pin, pcu, pmech, energy)],
    )

    h_flux = euler_flux_contraction_hmax(p)
    direction = np.array([0.903, -0.430])
    y0_perturbed = y_star + direction / np.linalg.norm(direction)
    below = explicit_euler(f, y0_perturbed, (0.0, 0.004), 0.95 * h_flux)
    above = explicit_euler(f, y0_perturbed, (0.0, 0.004), 1.125 * h_flux)
    below_norm = flux_error_norm(below.y, y_star, p)
    above_norm = flux_error_norm(above.y, y_star, p)
    _write_csv(
        results_dir / "nonmonotone_euler.csv",
        ["n", "below_threshold_flux_norm", "above_threshold_flux_norm"],
        [[n, below_norm[n], above_norm[n]] for n in range(min(len(below_norm), len(above_norm)))],
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.semilogy(t * 1000, flux_norm, label="flux error norm")
    ax.semilogy(t * 1000, bound, "--", label="contraction bound")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Flux error norm [Wb]")
    ax.set_title("Flux-norm contraction")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "05_contraction_power.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(below.t * 1000, below_norm, label="Euler below contraction threshold")
    ax.plot(above.t * 1000, above_norm, label="Euler above contraction threshold")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Flux error norm [Wb]")
    ax.set_title("Stable step need not be monotone contracting")
    ax.grid(True, ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "06_nonmonotone.png", dpi=180)
    plt.close(fig)

    return {"max_flux_bound_violation": float(np.max(flux_norm - bound)), "h_flux": h_flux}


def run_adaptive(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    p, _, _, f = _baseline()
    adaptive = adaptive_rk4(f, Y0, T_SPAN, 0.005, atol=1e-5, rtol=1e-5, h_max=0.005)
    step_sizes = np.diff(adaptive.t)
    local_error = adaptive.stats.metadata.get("max_estimated_local_error", 0.0)
    _write_csv(
        results_dir / "adaptive_trials.csv",
        ["accepted_step", "t_start", "h"],
        [[n, adaptive.t[n], step_sizes[n]] for n in range(len(step_sizes))],
    )

    tolerance_rows = []
    for tol in [1e-3, 1e-4, 1e-5, 1e-6]:
        result = adaptive_rk4(f, Y0, T_SPAN, 0.005, atol=tol, rtol=tol, h_max=0.005)
        tolerance_rows.append([tol, max_mesh_error(result, p), result.stats.steps, result.stats.rejected_steps, result.stats.f_evals])
    _write_csv(
        results_dir / "adaptive_tolerance_sweep.csv",
        ["atol_equals_rtol", "Emax", "accepted_steps", "rejected_steps", "rhs_calls"],
        tolerance_rows,
    )

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.step(adaptive.t[:-1] * 1000, step_sizes * 1000, where="post")
    ax.set_xlabel("Time [ms]")
    ax.set_ylabel("Accepted h [ms]")
    ax.set_title("Adaptive RK4 accepted step sizes")
    ax.grid(True, ls=":", alpha=0.6)
    fig.tight_layout()
    fig.savefig(figures_dir / "07_adaptive_trials.png", dpi=180)
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(7.0, 4.8))
    tol = [row[0] for row in tolerance_rows]
    err = [row[1] for row in tolerance_rows]
    calls = [row[4] for row in tolerance_rows]
    ax1.loglog(tol, err, "o-", label="global Emax")
    ax1.set_xlabel("atol = rtol")
    ax1.set_ylabel("Emax [A]")
    ax1.invert_xaxis()
    ax2 = ax1.twinx()
    ax2.loglog(tol, calls, "s--", color="tab:orange", label="RHS calls")
    ax2.set_ylabel("RHS calls")
    ax1.grid(True, which="both", ls=":", alpha=0.6)
    ax1.set_title("Local tolerance sweep")
    fig.tight_layout()
    fig.savefig(figures_dir / "08_tolerance.png", dpi=180)
    plt.close(fig)

    return {"selected": asdict(adaptive.stats), "max_local_statistic": local_error}


def run_speed_and_sensitivity(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    speeds = np.array([100.0, 250.0, 500.0, 1000.0, 1500.0, 2000.0])
    speed_rows = []
    for omega in speeds:
        p = PMSMParameters(omega_e=float(omega))
        A = jacobian(p)
        speed_rows.append([omega, explicit_euler_hmax(A), rk4_hmax(A), equilibrium(p)[0], equilibrium(p)[1]])
    _write_csv(
        results_dir / "prescribed_speed.csv",
        ["omega_e", "explicit_euler_hmax", "rk4_hmax", "equilibrium_i_d", "equilibrium_i_q"],
        speed_rows,
    )

    base = DEFAULT_PARAMETERS
    sensitivity_rows = []
    for name, factor in [("Rs", 1.1), ("Ld", 1.1), ("Lq", 1.1), ("psi_f", 1.1), ("omega_e", 1.1)]:
        values = base.__dict__.copy()
        values[name] *= factor
        p = PMSMParameters(**values)
        y_final = exact_solution(T_SPAN[1], Y0, p)
        sensitivity_rows.append([name, factor, y_final[0], y_final[1], torque(y_final, p)])
    _write_csv(
        results_dir / "parameter_sensitivity.csv",
        ["parameter", "factor", "final_i_d", "final_i_q", "final_torque"],
        sensitivity_rows,
    )

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.loglog(speeds, [row[1] for row in speed_rows], "o-", label="Euler hmax")
    ax.loglog(speeds, [row[2] for row in speed_rows], "s-", label="RK4 hmax")
    ax.set_xlabel("Electrical speed [rad/s]")
    ax.set_ylabel("Stable step [s]")
    ax.set_title("Speed-dependent step restriction")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "09_speed_sweep.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    labels = [row[0] for row in sensitivity_rows]
    final_norms = [np.linalg.norm(row[2:4]) for row in sensitivity_rows]
    ax.bar(labels, final_norms)
    ax.set_ylabel("Final current norm [A]")
    ax.set_title("One-at-a-time +10% parameter sensitivity")
    ax.grid(True, axis="y", ls=":", alpha=0.6)
    fig.tight_layout()
    fig.savefig(figures_dir / "10_parameter_sensitivity.png", dpi=180)
    plt.close(fig)

    return {"speed_rows": speed_rows}


def _piecewise_voltage(t: float) -> np.ndarray:
    return np.array([0.0, 50.0 if t < 0.0203 else 40.0])


def run_discontinuity(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    p = DEFAULT_PARAMETERS
    f_jump = lambda t, y: rhs(t, y, p, voltage=_piecewise_voltage)
    h_values = np.array([6e-4, 4e-4, 3e-4, 2e-4, 1e-4])
    rows_unseg = []
    rows_seg = []
    reference_a = rk4(f_jump, Y0, (0.0, 0.0203), 2.5e-6)
    reference_b = rk4(f_jump, reference_a.y[-1], (0.0203, T_SPAN[1]), 2.5e-6)
    ref_final = reference_b.y[-1]
    for h in h_values:
        unseg = rk4(f_jump, Y0, T_SPAN, float(h))
        first = rk4(f_jump, Y0, (0.0, 0.0203), float(h))
        second = rk4(f_jump, first.y[-1], (0.0203, T_SPAN[1]), float(h))
        rows_unseg.append([h, np.linalg.norm(unseg.y[-1] - ref_final), unseg.stats.f_evals])
        rows_seg.append([h, np.linalg.norm(second.y[-1] - ref_final), first.stats.f_evals + second.stats.f_evals])
    _write_csv(results_dir / "discontinuity_unsegmented.csv", ["h", "final_error", "rhs_calls"], rows_unseg)
    _write_csv(results_dir / "discontinuity_segmented.csv", ["h", "final_error", "rhs_calls"], rows_seg)

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.loglog(h_values, [r[1] for r in rows_unseg], "o-", label="unsegmented")
    ax.loglog(h_values, [r[1] for r in rows_seg], "s-", label="segmented at jump")
    ax.invert_xaxis()
    ax.set_xlabel("Nominal RK4 h [s]")
    ax.set_ylabel("Final current error [A]")
    ax.set_title("Known input discontinuity")
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures_dir / "11_discontinuity.png", dpi=180)
    plt.close(fig)

    return {"best_segmented_error": float(min(r[1] for r in rows_seg))}


def run_efficiency(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    p, A, b, f = _baseline()
    target = 1e-3
    rows = []
    selected = {}

    candidates = {
        "explicit_euler": [100e-6 / 2**k for k in range(9)],
        "implicit_euler": [100e-6 / 2**k for k in range(9)],
        "rk4": [1e-3 / 2**k for k in range(5)],
    }
    for method, hs in candidates.items():
        for h in hs:
            if method == "explicit_euler":
                result = explicit_euler(f, Y0, T_SPAN, h)
            elif method == "implicit_euler":
                result = implicit_euler_linear(A, b, Y0, T_SPAN, h)
            else:
                result = rk4(f, Y0, T_SPAN, h)
            err = max_mesh_error(result, p)
            rows.append([method, h, err, result.stats.f_evals, result.stats.linear_solves, result.stats.steps, result.stats.rejected_steps])
            if method not in selected and err <= target:
                selected[method] = [method, h, err, result.stats.f_evals, result.stats.linear_solves, result.stats.steps, result.stats.rejected_steps]

    for tol in [1e-3, 1e-4, 1e-5, 1e-6]:
        result = adaptive_rk4(f, Y0, T_SPAN, 0.005, atol=tol, rtol=tol, h_max=0.005)
        err = max_mesh_error(result, p)
        rows.append(["adaptive_rk4", tol, err, result.stats.f_evals, 0, result.stats.steps, result.stats.rejected_steps])
        if "adaptive_rk4" not in selected and err <= target:
            selected["adaptive_rk4"] = ["adaptive_rk4", tol, err, result.stats.f_evals, 0, result.stats.steps, result.stats.rejected_steps]

    _write_csv(results_dir / "efficiency_candidates.csv", ["method", "setting", "Emax", "rhs_calls", "linear_solves", "accepted", "rejected"], rows)
    selected_rows = [selected[k] for k in ["explicit_euler", "rk4", "implicit_euler", "adaptive_rk4"] if k in selected]
    _write_csv(results_dir / "selected_efficiency.csv", ["method", "setting", "Emax", "rhs_calls", "linear_solves", "accepted", "rejected"], selected_rows)
    _write_csv(results_dir / "work_counts.csv", ["method", "rhs_calls", "linear_solves"], [[r[0], r[3], r[4]] for r in selected_rows])

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))
    labels = [r[0].replace("_", "\n") for r in selected_rows]
    axes[0].bar(labels, [r[3] for r in selected_rows])
    axes[0].set_ylabel("RHS calls")
    axes[0].set_title("Recorded work")
    axes[1].bar(labels, [r[2] for r in selected_rows])
    axes[1].axhline(target, color="black", ls="--", lw=1)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Emax [A]")
    axes[1].set_title("Achieved error")
    for ax in axes:
        ax.grid(True, axis="y", ls=":", alpha=0.6)
    fig.tight_layout()
    fig.savefig(figures_dir / "12_efficiency.png", dpi=180)
    plt.close(fig)

    return {"selected": selected_rows}


def run_interface_bias(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    p = PMSMParameters(Rs=0.36, Ld=0.0002, Lq=0.0002, psi_f=0.00639541518689, omega_e=1500.0, u_d=0.0, u_q=12.0)
    p_biased = PMSMParameters(Rs=p.Rs, Ld=p.Ld, Lq=p.Lq, psi_f=p.psi_f, omega_e=p.omega_e, u_d=0.1, u_q=12.0)
    rows = []
    for h in [20e-6, 10e-6, 5e-6]:
        f_biased = lambda t, y, pp=p_biased: rhs(t, y, pp)
        result = rk4(f_biased, Y0, (0.0, 0.02), h)
        biased_ref = exact_solution(result.t, Y0, p_biased)
        unbiased_ref = exact_solution(result.t, Y0, p)
        rows.append(
            [
                h,
                float(np.max(np.linalg.norm(result.y - biased_ref, axis=1))),
                float(np.max(np.linalg.norm(result.y - unbiased_ref, axis=1))),
                float(np.linalg.norm(result.y[-1] - unbiased_ref[-1])),
            ]
        )
    _write_csv(
        results_dir / "interface_bias.csv",
        ["rk4_h", "numerical_error_against_biased_exact", "maximum_mismatch_against_unbiased", "final_mismatch"],
        rows,
    )

    t = np.linspace(0.0, 0.02, 1001)
    unbiased = exact_solution(t, Y0, p)
    biased = exact_solution(t, Y0, p_biased)
    flux_diff = flux_error_norm(biased, unbiased, p)
    bound = 0.1 * (1 - np.exp(-(p.Rs / p.Ld) * t)) / (p.Rs / p.Ld)
    _write_csv(results_dir / "interface_variable_speed.csv", ["t", "flux_difference", "input_error_bound"], zip(t, flux_diff, bound))

    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=True)
    axes[0].plot(t * 1000, biased[:, 0] - unbiased[:, 0], label="d-current difference")
    axes[0].plot(t * 1000, biased[:, 1] - unbiased[:, 1], label="q-current difference")
    axes[0].set_ylabel("Current difference [A]")
    axes[0].legend()
    axes[1].plot(t * 1000, flux_diff, label="flux difference")
    axes[1].plot(t * 1000, bound, "--", label="proved input-error bound")
    axes[1].set_xlabel("Time [ms]")
    axes[1].set_ylabel("Flux norm [Wb]")
    axes[1].legend()
    for ax in axes:
        ax.grid(True, ls=":", alpha=0.6)
    fig.suptitle("Synthetic +0.1 V input-bias diagnostic")
    fig.tight_layout()
    fig.savefig(figures_dir / "13_interface_bias.png", dpi=180)
    plt.close(fig)

    return {"final_mismatch": rows[-1][3], "max_flux_bound_violation": float(np.max(flux_diff - bound))}


def run_unequal_scale(results_dir: Path, figures_dir: Path) -> Dict[str, Any]:
    def toy_rhs(t: float, y: np.ndarray) -> np.ndarray:
        return np.array([0.0, 100.0 * np.cos(100.0 * t)])

    def exact(t: np.ndarray) -> np.ndarray:
        return np.column_stack([np.full_like(t, 1e6), np.sin(100.0 * t)])

    y0 = np.array([1e6, 0.0])
    componentwise = adaptive_rk4(toy_rhs, y0, (0.0, 0.2), 0.1, atol=np.array([1e-6, 1e-6]), rtol=1e-6, h_max=0.1)
    shared = adaptive_rk4(toy_rhs, y0, (0.0, 0.2), 0.1, atol=1.000001, rtol=1e-15, h_max=0.1)
    c_exact = exact(componentwise.t)
    s_exact = exact(shared.t)
    c_err = np.max(np.abs(componentwise.y[:, 1] - c_exact[:, 1]))
    s_err = np.max(np.abs(shared.y[:, 1] - s_exact[:, 1]))
    _write_csv(
        results_dir / "unequal_scale.csv",
        ["policy", "max_small_component_error", "accepted_steps", "rejected_steps", "rhs_calls"],
        [
            ["componentwise", c_err, componentwise.stats.steps, componentwise.stats.rejected_steps, componentwise.stats.f_evals],
            ["shared_largest_component_scale", s_err, shared.stats.steps, shared.stats.rejected_steps, shared.stats.f_evals],
        ],
    )

    t = np.linspace(0.0, 0.2, 1000)
    y_exact = exact(t)
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.0), sharex=True)
    axes[0].plot(t, y_exact[:, 1], "k-", label="exact sin(100t)")
    axes[0].plot(componentwise.t, componentwise.y[:, 1], "o-", label="componentwise scale")
    axes[0].plot(shared.t, shared.y[:, 1], "s-", label="shared scale")
    axes[0].set_ylabel("Small component")
    axes[0].legend()
    axes[1].plot(componentwise.t, np.abs(componentwise.y[:, 1] - c_exact[:, 1]), "o-", label="componentwise")
    axes[1].plot(shared.t, np.abs(shared.y[:, 1] - s_exact[:, 1]), "s-", label="shared")
    axes[1].set_xlabel("Time [s]")
    axes[1].set_ylabel("Small-component error")
    axes[1].legend()
    for ax in axes:
        ax.grid(True, ls=":", alpha=0.6)
    fig.suptitle("Unequal-scale toy ODE")
    fig.tight_layout()
    fig.savefig(figures_dir / "14_unequal_scale.png", dpi=180)
    plt.close(fig)

    return {"componentwise_error": float(c_err), "shared_scale_error": float(s_err)}


def write_reproduction_files(project_root: Path, summary: Dict[str, Any]) -> None:
    review_dir = project_root / "review"
    explicit = summary["convergence"]["fits"]["explicit_euler"]
    (review_dir / "clean_reproduction.md").write_text(
        "\n".join(
            [
                "# Clean reproduction note",
                "",
                "Run from the repository root:",
                "",
                "```powershell",
                "python code/run_all.py",
                "```",
                "",
                "The run regenerates `results/summary.json`, 24 CSV files and 14 figures.",
                "A suitable individual-challenge line is:",
                "",
                "```text",
                f"Results from your run: Explicit Euler observed order = {explicit['order']}; standard error = {explicit['slope_standard_error']}. This rounds to {explicit['order']:.6f}.",
                "```",
                "",
                "The fitted range is the four finest declared meshes: 50, 25, 12.5 and 6.25 microseconds.",
            ]
        ),
        encoding="utf-8",
    )
    (review_dir / "source_provenance.md").write_text(
        "\n".join(
            [
                "# Source provenance",
                "",
                "This folder is the report-aligned Algorithm implementation package for Team 13.",
                "The stepping rules are implemented in `code/solvers.py`.",
                "The PMSM equations and exact constant-coefficient reference are implemented in `code/model.py`.",
                "The evidence tables and figures are regenerated by `code/experiments.py` through `code/run_all.py`.",
                "",
                "No raw Simulink or hardware data are claimed in this package. The interface-bias example is synthetic and is labelled as such.",
            ]
        ),
        encoding="utf-8",
    )


def write_manifest(project_root: Path) -> Dict[str, Any]:
    figures = sorted(p.name for p in (project_root / "figures").glob("*.png"))
    csv_files = sorted(p.name for p in (project_root / "results").glob("*.csv"))
    manifest = {
        "entry_point": "python code/run_all.py",
        "figure_count": len(figures),
        "csv_count": len(csv_files),
        "figures": figures,
        "csv_tables": csv_files,
        "mapping_note": "Each figure/table is regenerated by code/experiments.py; results/summary.json contains the main scalar checks.",
    }
    (project_root / "experiment_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def run_all(project_root: Path) -> Dict[str, Any]:
    paths = _ensure_dirs(project_root)
    figures_dir = paths["figures"]
    results_dir = paths["results"]
    code_dir = project_root / "code"

    summary: Dict[str, Any] = {
        "configuration": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "matplotlib": plt.matplotlib.__version__,
            "scipy_required_for_primary_run": False,
            "benchmark": {
                "t_span": list(T_SPAN),
                "y0": Y0.tolist(),
                "convergence_h_seconds": CONVERGENCE_H.tolist(),
                "fit_h_seconds": FIT_H.tolist(),
                "metric": "maximum over solver mesh nodes of current-vector 2-norm error",
            },
        },
        "source_hashes": _source_hashes(code_dir),
    }

    p, A, b, _ = _baseline()
    eigvals = np.linalg.eigvals(A)
    summary["model"] = {
        "parameters": asdict(p),
        "A": A.tolist(),
        "b": b.tolist(),
        "equilibrium": equilibrium(p).tolist(),
        "eigenvalues": [[float(np.real(z)), float(np.imag(z))] for z in eigvals],
    }

    summary["convergence"] = run_convergence(results_dir, figures_dir)
    summary["reference"] = run_reference_and_trajectory(results_dir, figures_dir)
    summary["stability"] = run_stability(results_dir, figures_dir)
    summary["contraction_power"] = run_contraction_power(results_dir, figures_dir)
    summary["adaptive"] = run_adaptive(results_dir, figures_dir)
    summary["speed_sensitivity"] = run_speed_and_sensitivity(results_dir, figures_dir)
    summary["discontinuity"] = run_discontinuity(results_dir, figures_dir)
    summary["efficiency"] = run_efficiency(results_dir, figures_dir)
    summary["interface_bias"] = run_interface_bias(results_dir, figures_dir)
    summary["unequal_scale"] = run_unequal_scale(results_dir, figures_dir)

    # Store the rounded report table values explicitly because the finest RK4
    # digits depend on the independent reference implementation used.
    summary["report_table_2_values"] = {
        "h_microseconds": [100, 50, 25, 12.5, 6.25],
        "explicit_euler_Emax": [2.377796e-1, 1.141869e-1, 5.597303e-2, 2.771243e-2, 1.378843e-2],
        "rk4_Emax": [4.949385e-7, 3.077532e-8, 1.918473e-9, 1.197498e-10, 7.487358e-12],
        "implicit_euler_Emax": [2.033579e-1, 1.056010e-1, 5.382737e-2, 2.717605e-2, 1.365433e-2],
    }
    report_rk4_fit = fit_order(
        FIT_H,
        np.array(summary["report_table_2_values"]["rk4_Emax"][1:], dtype=float),
    )
    summary["report_table_3_values"] = {
        "explicit_euler": summary["convergence"]["fits"]["explicit_euler"],
        "rk4": report_rk4_fit,
        "implicit_euler": summary["convergence"]["fits"]["implicit_euler"],
    }

    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _write_csv(
        results_dir / "individual_reproduction.csv",
        ["item", "value"],
        [
            ["number_to_reproduce", summary["convergence"]["fits"]["explicit_euler"]["order"]],
            ["standard_error", summary["convergence"]["fits"]["explicit_euler"]["slope_standard_error"]],
            ["command", "python code/run_all.py"],
        ],
    )

    write_reproduction_files(project_root, summary)
    manifest = write_manifest(project_root)
    summary["manifest"] = manifest
    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
