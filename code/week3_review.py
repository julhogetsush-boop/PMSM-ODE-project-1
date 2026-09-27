"""Tutorial 3 quantitative review: one observation -> one defensible figure.

Chosen finding: the Week 2 convergence figure reported fitted slopes but did
not show the theoretical slope, fit uncertainty, or the round-off regime.
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np

from methods import solve_fixed
from model import exact_solution, load_parameters, make_model

ROOT = Path(__file__).resolve().parents[1]


def run_week3(figures=None, results=None):
    figures = ROOT / "figures" if figures is None else Path(figures)
    results = ROOT / "results" if results is None else Path(results)
    figures.mkdir(exist_ok=True)
    results.mkdir(exist_ok=True)

    parameters = load_parameters()
    rhs, jacobian = make_model(parameters)
    final_time = parameters["T"]
    reference = exact_solution(final_time, parameters["y0"], parameters)

    # Powers of two show both the order-4 range and the floating-point floor.
    step_counts = np.array(
        [60, 120, 240, 480, 960, 1920, 3840, 7680,
         15360, 30720, 61440, 122880],
        dtype=int,
    )
    step_sizes = final_time / step_counts
    errors = []
    for count, step in zip(step_counts, step_sizes):
        _, states, _ = solve_fixed(
            rhs, parameters["y0"], 0.0, final_time, step,
            method="rk4", jac=jacobian,
        )
        errors.append(float(np.max(np.abs(states[-1] - reference))))
    errors = np.array(errors)

    # The coarsest point is pre-asymptotic; the last four approach round-off.
    fit_mask = (step_counts >= 120) & (step_counts <= 7680)
    log_h = np.log(step_sizes[fit_mask])
    log_error = np.log(errors[fit_mask])
    coefficients, covariance = np.polyfit(log_h, log_error, 1, cov=True)
    slope, intercept = coefficients
    slope_standard_error = float(np.sqrt(covariance[0, 0]))

    fit_h = np.logspace(np.log10(step_sizes.min()), np.log10(step_sizes.max()), 300)
    fit_log_h = np.log(fit_h)
    fit_log_error = slope * fit_log_h + intercept
    fit_standard_error = np.sqrt(
        covariance[1, 1]
        + 2 * fit_log_h * covariance[0, 1]
        + fit_log_h**2 * covariance[0, 0]
    )

    # Anchor the order-4 benchmark at the geometric centre of the fit range.
    anchor_h = float(np.exp(np.mean(log_h)))
    anchor_error = float(np.exp(slope * np.log(anchor_h) + intercept))
    theoretical_error = anchor_error * (fit_h / anchor_h) ** 4
    floor_boundary = final_time / 30720

    rows = [
        {
            "N": int(count),
            "h_s": float(step),
            "endpoint_Linf_error_A": float(error),
            "used_in_fit": bool(use),
            "roundoff_region": bool(step <= floor_boundary),
        }
        for count, step, error, use in zip(step_counts, step_sizes, errors, fit_mask)
    ]
    with (results / "week3_rk4_convergence.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.loglog(
        step_sizes, errors, "o-", color="#0072B2", markerfacecolor="white",
        label="RK4 endpoint error",
    )
    ax.loglog(
        fit_h, theoretical_error, "k--", linewidth=1.4,
        label="theoretical order 4",
    )
    ax.loglog(
        fit_h, np.exp(fit_log_error), color="#D55E00", linestyle="-.",
        label=f"fit: {slope:.3f} +/- {slope_standard_error:.3f}",
    )
    ax.fill_between(
        fit_h,
        np.exp(fit_log_error - 1.96 * fit_standard_error),
        np.exp(fit_log_error + 1.96 * fit_standard_error),
        color="#D55E00", alpha=0.16, label="95% regression band",
    )
    ax.axvspan(step_sizes.min(), floor_boundary, color="grey", alpha=0.16)
    ax.axvline(floor_boundary, color="grey", linestyle=":", linewidth=1)
    ax.annotate(
        "round-off-dominated\n(non-monotone error)",
        xy=(floor_boundary, errors[9]),
        xytext=(floor_boundary * 6, errors[9] * 7),
        arrowprops={"arrowstyle": "->", "color": "0.35"},
        fontsize=9,
    )
    ax.set_xlabel("step size h [s]")
    ax.set_ylabel("endpoint current error, L-infinity [A]")
    ax.set_title("RK4 achieves order 4 before the floating-point floor")
    ax.grid(True, which="both", linestyle=":", alpha=0.35)
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figures / "week3_rk4_quantitative.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "feedback_item": "Convergence plot missing theoretical slope, uncertainty, and round-off annotation",
        "method": "RK4",
        "error_metric": "endpoint componentwise L-infinity current error [A]",
        "fit_N_min": 120,
        "fit_N_max": 7680,
        "fit_h_max_s": float(final_time / 120),
        "fit_h_min_s": float(final_time / 7680),
        "slope": float(slope),
        "slope_standard_error": slope_standard_error,
        "roundoff_boundary_h_s": float(floor_boundary),
        "smallest_reported_error_A": float(errors.min()),
        "band_interpretation": "deterministic regression diagnostic, not repeated-sample uncertainty",
    }
    (results / "week3_review.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # A narrow, meaningful check: this revision must support the stated claim.
    assert abs(slope - 4) < 0.1
    assert slope_standard_error < 0.05
    assert errors[10] >= errors[9]  # the selected floor is visibly non-monotone
    return summary


if __name__ == "__main__":
    print(json.dumps(run_week3(), indent=2, ensure_ascii=False))
