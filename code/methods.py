"""Model-agnostic Euler, RK4, Newton and adaptive implicit Euler.

All states are 1-D NumPy arrays. Rejected trials never advance time.
"""
import numpy as np


class NewtonFailure(RuntimeError):
    pass


def newton(F, JF, w0, tol=1e-12, max_iter=40, damped=False, norm=np.inf):
    """Return root and history; k counts accepted Newton updates, not probes."""
    if not np.isfinite(tol) or tol <= 0 or max_iter < 1:
        raise ValueError("Positive tolerance and iteration limit required")
    w = np.asarray(w0, dtype=float).copy()
    history = []
    for k in range(max_iter+1):
        with np.errstate(over="ignore", invalid="ignore"):
            residual = np.asarray(F(w))
            r = float(np.linalg.norm(residual, ord=norm))
        history.append({"iteration": k, "residual": r, "alpha": 0.0})
        if not np.isfinite(r):
            raise NewtonFailure("Nonfinite Newton residual")
        if r < tol:
            return w, history
        if k == max_iter:
            break
        try:
            delta = np.linalg.solve(JF(w), -residual)
        except np.linalg.LinAlgError as exc:
            raise NewtonFailure("Singular Newton Jacobian") from exc
        alpha = 1.0
        if damped:
            for _ in range(40):
                with np.errstate(over="ignore", invalid="ignore"):
                    trial_r = np.linalg.norm(F(w+alpha*delta), ord=norm)
                if np.isfinite(trial_r) and trial_r < r:
                    break
                alpha *= 0.5
            else:
                raise NewtonFailure("Backtracking failed to decrease residual")
        w = w + alpha*delta
        history[-1]["alpha"] = alpha
    raise NewtonFailure(f"Newton exceeded {max_iter} updates")


def euler_step(f, t, y, h):
    return y + h*f(t, y)


def rk4_step(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t+h/2, y+h*k1/2)
    k3 = f(t+h/2, y+h*k2/2)
    k4 = f(t+h, y+h*k3)
    return y + h*(k1+2*k2+2*k3+k4)/6


def implicit_step(f, jac, t, y, h, tol=1e-12, damped=False, max_iter=40):
    F = lambda w: w-y-h*f(t+h, w)
    JF = lambda w: np.eye(y.size)-h*jac(t+h, w)
    return newton(F, JF, y, tol=tol, damped=damped, max_iter=max_iter)


def time_grid(t0, T, h):
    if not np.all(np.isfinite([t0, T, h])) or h <= 0 or T <= t0:
        raise ValueError("Require finite t0 < T and h > 0")
    # Construct from n*h, avoiding cumulative floating point time drift.
    n = int(np.ceil((T-t0)/h))
    if n > 2_000_000:
        raise ValueError("Too many steps; review h and T")
    t = t0+h*np.arange(n+1, dtype=float)
    t = t[t < T-8*np.finfo(float).eps*max(1, abs(T))]
    return np.r_[t, T]


def solve_fixed(f, y0, t0, T, h, method="euler", jac=None, **newton_options):
    t = time_grid(t0, T, h)
    y0 = np.atleast_1d(np.asarray(y0, dtype=float))
    if y0.ndim != 1 or not np.all(np.isfinite(y0)):
        raise ValueError("Initial state must be a finite vector")
    if method not in ("euler", "rk4", "implicit"):
        raise ValueError(method)
    if method == "implicit" and jac is None:
        raise ValueError("Implicit Euler requires a Jacobian")
    Y = np.empty((len(t), len(y0)))
    Y[0] = y0
    updates = []
    for n, step in enumerate(np.diff(t)):
        with np.errstate(over="ignore", invalid="ignore"):
            if method == "implicit":
                Y[n+1], hist = implicit_step(f, jac, t[n], Y[n], step, **newton_options)
                updates.append(len(hist)-1)
            else:
                func = euler_step if method == "euler" else rk4_step
                Y[n+1] = func(f, t[n], Y[n], step)
        if not np.all(np.isfinite(Y[n+1])):
            raise FloatingPointError(f"Nonfinite state at t={t[n+1]:.9g}")
    return t, Y, {"steps": len(t)-1, "newton_updates": sum(updates),
                  "max_newton_updates": max(updates, default=0)}


def adaptive_implicit(f, jac, y0, t0, T, h0, atol=1e-5, rtol=1e-4,
                      hmax=None, hmin=1e-12, max_trials=100000, **options):
    """Step doubling: accept two half steps; p=1, error denominator 2**p-1=1.

    Newton uses state residual units; local control uses component scales.
    All failed/rejected trials are logged. No Richardson correction is used.
    """
    time_grid(t0, T, h0)  # validate interval and starting step
    y = np.atleast_1d(np.asarray(y0, dtype=float)).copy()
    atol = np.broadcast_to(np.asarray(atol, dtype=float), y.shape)
    hmax = T-t0 if hmax is None else hmax
    if (not np.all(np.isfinite(y)) or not np.all(np.isfinite(atol)) or np.any(atol <= 0)
            or not np.isfinite(rtol) or rtol < 0 or not 0 < hmin <= hmax):
        raise ValueError("Invalid state, tolerances or step bounds")
    t, h = float(t0), min(h0, hmax)
    times, states, trials = [t], [y.copy()], []
    for _ in range(max_trials):
        if t >= T:
            return np.array(times), np.array(states), trials
        h = min(h, hmax, T-t)
        if h < hmin or t+h == t:
            raise RuntimeError("Adaptive step fell below hmin")
        error, count, reason = np.inf, 0, ""
        try:
            coarse, hist = implicit_step(f, jac, t, y, h, **options)
            count += len(hist)-1
            half, hist = implicit_step(f, jac, t, y, h/2, **options)
            count += len(hist)-1
            fine, hist = implicit_step(f, jac, t+h/2, half, h/2, **options)
            count += len(hist)-1
            scale = atol + rtol*np.maximum(abs(y), abs(fine))
            error = float(np.max(abs(fine-coarse)/scale))
        except (NewtonFailure, FloatingPointError) as exc:
            reason = str(exc)
        accepted = bool(np.isfinite(error) and error <= 1)
        trials.append({"t": t, "h": h, "error_ratio": error, "accepted": accepted,
                       "newton_updates_completed_solves": count, "failure": reason})
        if accepted:
            t = min(T, t+h)
            y = fine
            times.append(t)
            states.append(y.copy())
        # A rejection leaves t and y untouched.
        factor = 2.0 if error == 0 else np.clip(0.9*error**(-0.5), 0.2, 2.0)
        h *= factor if accepted else min(0.5, factor)
    raise RuntimeError("Adaptive trial limit exceeded")


def finite_difference_jacobian(f, t, y):
    y = np.asarray(y, dtype=float)
    J = np.empty((len(y), len(y)))
    for j in range(len(y)):
        eps = np.cbrt(np.finfo(float).eps)*max(1.0, abs(y[j]))
        e = np.zeros_like(y)
        e[j] = eps
        J[:, j] = (f(t, y+e)-f(t, y-e))/(2*eps)
    return J
