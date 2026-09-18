"""前两周课堂例题：Logistic / Van der Pol / Robertson / RC diode.

Equations and experiment settings follow the supplied course materials.
Results are computed here, never copied from the slides.
"""
import csv
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from methods import newton, solve_fixed, euler_step


def save_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def reference(f, y0, T, method="Radau"):
    """Independent solver, checked by tightening its tolerances."""
    base = solve_ivp(f, (0, T), y0, method=method, rtol=1e-10,
                     atol=1e-12, dense_output=True)
    tight = solve_ivp(f, (0, T), y0, method=method, rtol=1e-12,
                      atol=1e-14, dense_output=True)
    if not (base.success and tight.success):
        raise RuntimeError("Reference solve failed")
    tt = np.linspace(0, T, 2001)
    gap = float(np.max(abs(base.sol(tt)-tight.sol(tt))))
    if gap > 1e-7:
        raise RuntimeError(f"Reference not resolved: {gap}")
    return tight, gap


def logistic(t, y):
    return y*(1-y/10)


def vanderpol(mu):
    return lambda t, y: np.array([y[1], mu*(1-y[0]**2)*y[1]-y[0]])


def robertson(t, y):
    a, b, c = y
    return np.array([-0.04*a+1e4*b*c, 0.04*a-1e4*b*c-3e7*b*b, 3e7*b*b])


def robertson_jac(t, y):
    a, b, c = y
    return np.array([[-0.04, 1e4*c, 1e4*b],
                     [0.04, -1e4*c-6e7*b, -1e4*b], [0, 6e7*b, 0]])


def rc_residual(v, old, h):
    # Current residual in A, as in Tutorial 2; not the voltage residual.
    v1, v2 = v
    return np.array([1e-5*(v1-old[0])/h-(5-v1)/100+(v1-v2)/100,
                     1e-5*(v2-old[1])/h-(v1-v2)/100+1e-12*np.expm1(v2/0.02585)])


def rc_jac(v, h):
    return np.array([[1e-5/h+0.02, -0.01],
                     [-0.01, 1e-5/h+0.01+1e-12/0.02585*np.exp(v[1]/0.02585)]])


def rc_step(old, h, damped):
    return newton(lambda v: rc_residual(v, old, h), lambda v: rc_jac(v, h),
                  old, tol=1e-12, max_iter=300, damped=damped, norm=2)


def run_labs(figures, results):
    summary = {}
    # Week 1: two hand-computed steps, then a genuine convergence experiment.
    one = euler_step(logistic, 0, np.array([0.5]), 0.2)
    two = euler_step(logistic, 0.2, one, 0.2)
    # Slide typo: 0.595*(1-0.595/10)=0.5595975, hence y2=0.7069195.
    assert np.allclose([one[0], two[0]], [0.595, 0.7069195], atol=1e-14, rtol=0)
    hs = np.array([0.4, 0.2, 0.1, 0.05, 0.025])
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4))
    for h in hs:
        t, Y, _ = solve_fixed(logistic, [0.5], 0, 10, h)
        exact = 10/(1+19*np.exp(-t))
        rows.append({"h": h, "endpoint_error": abs(Y[-1, 0]-exact[-1]),
                     "max_grid_error": float(np.max(abs(Y[:, 0]-exact)))})
    axes[0].plot(t, Y[:, 0], label="Euler, h=0.025")
    axes[0].plot(t, exact, "k--", label="Exact")
    axes[0].set(xlabel="t", ylabel="y", title="Logistic calibration")
    axes[0].legend()
    errors = np.array([r["endpoint_error"] for r in rows])
    order = float(np.polyfit(np.log(hs[-3:]), np.log(errors[-3:]), 1)[0])
    axes[1].loglog(hs, errors, "o-", label=f"Euler, last-3 slope={order:.3f}")
    axes[1].set_xticks(sorted(hs), labels=[f"{v:g}" for v in sorted(hs)])
    axes[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axes[1].set(xlabel="h", ylabel="Endpoint absolute error", title="Same T=10 for every run")
    axes[1].legend()
    save_figure(fig, figures/"lab_logistic.png")
    save_csv(results/"lab_logistic.csv", rows)
    summary["logistic_order_last3"] = order
    # Scalar stability: h=2 is bounded but does not decay.
    rows = []
    for h in [1.9, 2.0, 2.1]:
        t, Y, _ = solve_fixed(lambda t, y: -y, [1], 0, 200, h)
        rows.append({"h": h, "amplification_magnitude": abs(1-h),
                     "abs_y_T": abs(Y[-1, 0])})
    save_csv(results/"lab_scalar_stability.csv", rows)
    # Tutorial 1: use the prescribed initial data, not the failure-case data.
    f = vanderpol(1)
    oracle, gap = reference(f, [0.5, 0], 10, "DOP853")
    rows = []
    for h in [0.1, 0.05, 0.025, 0.0125, 0.00625]:
        t, Y, _ = solve_fixed(f, [0.5, 0], 0, 10, h)
        rows.append({"h": h, "endpoint_L2_error": float(np.linalg.norm(Y[-1]-oracle.sol(10))),
                     "max_grid_Linf_error": float(np.max(abs(Y-oracle.sol(t).T)))})
    save_csv(results/"lab_vanderpol_week1.csv", rows)
    summary["vdp_week1_y1_T"] = float(oracle.sol(10)[0])
    summary["vdp_week1_reference_gap"] = gap
    summary["vdp_week1_order_prescribed4"] = float(np.polyfit(
        np.log([r["h"] for r in rows[:4]]), np.log([r["endpoint_L2_error"] for r in rows[:4]]), 1)[0])
    # Week 2: reproduce a failed explicit run, while recording finite output.
    f_bad = vanderpol(10)
    ref_bad, gap = reference(f_bad, [2, 0], 20)
    tt = np.linspace(0, 20, 2001)
    ts, ys = [0.0], [np.array([2.0, 0.0])]
    for n in range(200):
        trial = euler_step(f_bad, ts[-1], ys[-1], 0.1)
        ts.append((n+1)*0.1)
        ys.append(trial)
        if np.max(abs(trial)) > 1e6:
            break
    Ybad = np.array(ys)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4))
    phase = oracle.sol(np.linspace(0, 10, 1500))
    axes[0].plot(phase[0], phase[1])
    axes[0].set(xlabel="y1", ylabel="y2", title="Tutorial 1: mu=1 reference phase plane")
    axes[1].semilogy(tt, np.maximum(abs(ref_bad.sol(tt)[0]), 1e-12), label="Radau reference")
    axes[1].semilogy(ts, np.maximum(abs(Ybad[:, 0]), 1e-12), "o-", ms=2, label="Euler h=0.1; cutoff at 1e6")
    axes[1].set(xlabel="t", ylabel="|y1|", title="Week 2: mu=10 explicit instability")
    axes[1].legend(fontsize=8)
    save_figure(fig, figures/"lab_vanderpol.png")
    summary["vdp_failure_last_t"] = ts[-1]
    summary["vdp_failure_peak_state"] = float(np.max(abs(Ybad)))
    summary["vdp_failure_reference_gap"] = gap
    # Robertson: retain the full coupled Jacobian and measure mass separately.
    robref, gap = reference(robertson, [1, 0, 0], 40)
    t, Y, counts = solve_fixed(robertson, [1, 0, 0], 0, 40, 0.1,
                               "implicit", robertson_jac, max_iter=300)
    errors = abs(Y[-1]-robref.sol(40))
    mass_error = float(np.max(abs(Y.sum(axis=1)-1)))
    save_csv(results/"lab_robertson.csv", [{"component": j+1, "implicit_T": Y[-1,j],
              "reference_T": robref.sol(40)[j], "absolute_error": errors[j]} for j in range(3)])
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4))
    for j in range(3):
        axes[0].plot(t, Y[:,j], label=f"y{j+1}, IE h=0.1")
        axes[0].plot(t, robref.sol(t)[j], "k--", lw=0.7)
    axes[0].set(xlabel="t", ylabel="Concentration", title="Robertson: dashed curves are Radau")
    axes[0].legend(fontsize=8)
    axes[1].semilogy(t, np.maximum(abs(Y.sum(axis=1)-1), 1e-18))
    axes[1].set(xlabel="t", ylabel="|sum(y)-1| (display floor 1e-18)", title="Mass defect is not trajectory error")
    save_figure(fig, figures/"lab_robertson.png")
    summary["robertson_mass_error"] = mass_error
    summary["robertson_endpoint_error_max"] = float(max(errors))
    summary["robertson_reference_gap"] = gap
    # RC challenge: current residual, 2-norm, old-state guess, exact endpoint.
    fig, ax = plt.subplots(figsize=(7, 4))
    for damped in [False, True]:
        label = "damped" if damped else "full"
        root, hist = rc_step(np.zeros(2), 1e-3, damped)
        ax.semilogy([r["iteration"] for r in hist], [r["residual"] for r in hist], "o-", label=label)
        save_csv(results/f"lab_rc_trace_{label}.csv", hist)
        summary[f"rc_first_step_{label}_updates"] = len(hist)-1
        summary[f"rc_first_step_{label}_root"] = root.tolist()
    ax.set(xlabel="Accepted Newton update k", ylabel="Current residual 2-norm (A)",
           title="RC diode: same h=1 ms, different nonlinear iterations")
    ax.legend()
    save_figure(fig, figures/"lab_rc_newton.png")
    rows = []
    for h_ms in [0.1, 0.2, 0.4, 0.8, 1, 2, 4, 8, 20]:
        row = {"h_ms": h_ms}
        for damped in [False, True]:
            label = "damped" if damped else "full"
            old, t, maxit = np.zeros(2), 0.0, 0
            while t < 0.02-1e-15:
                step = min(h_ms*1e-3, 0.02-t)
                old, hist = rc_step(old, step, damped)
                maxit = max(maxit, len(hist)-1)
                t += step
            row[f"{label}_max_updates"] = maxit
            row[f"{label}_passes_10"] = maxit <= 10
        rows.append(row)
    save_csv(results/"lab_rc_step_sweep.csv", rows)
    for label in ["full", "damped"]:
        summary[f"rc_{label}_hmax_ms"] = max(r["h_ms"] for r in rows if r[f"{label}_passes_10"])
    assert summary["rc_first_step_full_updates"] == 21
    assert summary["rc_first_step_damped_updates"] == 6
    assert mass_error < 1e-9
    return summary
