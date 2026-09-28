"""Numerical time integrators used in the PMSM ODE report.

The implementations are intentionally small and explicit. They are the
student-facing algorithms, while NumPy is only used for array arithmetic and
linear solves.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray


RHS = Callable[[float, NDArray[np.float64]], NDArray[np.float64]]
Jacobian = Callable[[float, NDArray[np.float64]], NDArray[np.float64]]


@dataclass(frozen=True)
class SolverStats:
    steps: int
    f_evals: int
    rejected_steps: int = 0
    linear_solves: int = 0
    newton_iterations: int = 0
    metadata: Dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class SolverResult:
    t: NDArray[np.float64]
    y: NDArray[np.float64]
    method: str
    stats: SolverStats


def _validate_span(t_span: Tuple[float, float], h: float) -> Tuple[float, float]:
    t0, t1 = map(float, t_span)
    if not np.isfinite([t0, t1, h]).all():
        raise ValueError("t0, t1 and h must be finite.")
    if t1 < t0:
        raise ValueError("t1 must be at least t0.")
    if h <= 0:
        raise ValueError("h must be positive.")
    return t0, t1


def time_grid(t_span: Tuple[float, float], h: float) -> NDArray[np.float64]:
    """Return a clipped fixed-step grid including the final time."""
    t0, t1 = _validate_span(t_span, h)
    n_steps = int(np.ceil((t1 - t0) / h))
    t = np.empty(n_steps + 1, dtype=float)
    t[0] = t0
    for n in range(n_steps):
        t[n + 1] = min(t[n] + h, t1)
    return t


def explicit_euler(f: RHS, y0: ArrayLike, t_span: Tuple[float, float], h: float) -> SolverResult:
    """Solve y' = f(t, y) using explicit Euler."""
    t = time_grid(t_span, h)
    y = np.empty((len(t),) + np.shape(y0), dtype=float)
    y[0] = np.asarray(y0, dtype=float)

    for n in range(len(t) - 1):
        step = t[n + 1] - t[n]
        y[n + 1] = y[n] + step * f(t[n], y[n])
        if not np.isfinite(y[n + 1]).all():
            raise FloatingPointError("explicit Euler produced a non-finite state.")

    steps = len(t) - 1
    return SolverResult(t, y, "explicit_euler", SolverStats(steps, steps))


def rk4_step(f: RHS, t: float, y: NDArray[np.float64], h: float) -> NDArray[np.float64]:
    """One classical fourth-order Runge-Kutta step."""
    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2)
    k4 = f(t + h, y + h * k3)
    y_next = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    if not np.isfinite(y_next).all():
        raise FloatingPointError("RK4 produced a non-finite state.")
    return y_next


def rk4(f: RHS, y0: ArrayLike, t_span: Tuple[float, float], h: float) -> SolverResult:
    """Solve y' = f(t, y) using fixed-step classical RK4."""
    t = time_grid(t_span, h)
    y = np.empty((len(t),) + np.shape(y0), dtype=float)
    y[0] = np.asarray(y0, dtype=float)

    for n in range(len(t) - 1):
        y[n + 1] = rk4_step(f, t[n], y[n], t[n + 1] - t[n])

    steps = len(t) - 1
    return SolverResult(t, y, "rk4", SolverStats(steps, 4 * steps))


def implicit_euler_newton(
    f: RHS,
    jac: Jacobian,
    y0: ArrayLike,
    t_span: Tuple[float, float],
    h: float,
    *,
    residual_atol: float = 1e-11,
    residual_rtol: float = 1e-11,
    max_iter: int = 20,
) -> SolverResult:
    """Solve y' = f(t, y) with implicit Euler and Newton corrections."""
    t = time_grid(t_span, h)
    y = np.empty((len(t),) + np.shape(y0), dtype=float)
    y[0] = np.asarray(y0, dtype=float)
    f_evals = 0
    linear_solves = 0
    newton_iterations = 0

    for n in range(len(t) - 1):
        step = t[n + 1] - t[n]
        t_next = t[n + 1]
        w = y[n].copy()
        scale = residual_atol + residual_rtol * max(np.linalg.norm(w, np.inf), np.linalg.norm(y[n], np.inf))

        for _ in range(max_iter):
            F = w - y[n] - step * f(t_next, w)
            f_evals += 1
            residual = np.linalg.norm(F, ord=np.inf)
            if residual <= scale:
                break

            matrix = np.eye(w.size) - step * jac(t_next, w)
            delta = np.linalg.solve(matrix, -F)
            linear_solves += 1
            newton_iterations += 1

            trial_factor = 1.0
            accepted_trial = False
            while trial_factor >= 2.0 ** -12:
                candidate = w + trial_factor * delta
                F_candidate = candidate - y[n] - step * f(t_next, candidate)
                f_evals += 1
                if np.linalg.norm(F_candidate, ord=np.inf) <= (1.0 - 1e-4 * trial_factor) * residual:
                    w = candidate
                    accepted_trial = True
                    break
                trial_factor *= 0.5
            if not accepted_trial:
                raise RuntimeError("Newton backtracking failed to reduce the residual.")
        else:
            raise RuntimeError("implicit Euler Newton iteration failed to converge.")

        if not np.isfinite(w).all():
            raise FloatingPointError("implicit Euler produced a non-finite state.")
        y[n + 1] = w

    steps = len(t) - 1
    return SolverResult(
        t,
        y,
        "implicit_euler_newton",
        SolverStats(
            steps,
            f_evals,
            linear_solves=linear_solves,
            newton_iterations=newton_iterations,
        ),
    )


def implicit_euler_linear(
    A: ArrayLike,
    b: ArrayLike,
    y0: ArrayLike,
    t_span: Tuple[float, float],
    h: float,
) -> SolverResult:
    """Specialized implicit Euler for y' = A y + b."""
    t = time_grid(t_span, h)
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    y = np.empty((len(t),) + np.shape(y0), dtype=float)
    y[0] = np.asarray(y0, dtype=float)

    for n in range(len(t) - 1):
        step = t[n + 1] - t[n]
        y[n + 1] = np.linalg.solve(np.eye(A.shape[0]) - step * A, y[n] + step * b)

    steps = len(t) - 1
    return SolverResult(
        t,
        y,
        "implicit_euler_linear",
        SolverStats(steps, 0, linear_solves=steps, newton_iterations=steps),
    )


def adaptive_rk4(
    f: RHS,
    y0: ArrayLike,
    t_span: Tuple[float, float],
    h_initial: float,
    *,
    rtol: float = 1e-6,
    atol: float | ArrayLike = 1e-9,
    h_min: float = 1e-14,
    h_max: float | None = None,
    safety: float = 0.9,
    min_factor: float = 0.2,
    max_factor: float = 5.0,
    max_trials: int = 200000,
) -> SolverResult:
    """Adaptive RK4 using step doubling and componentwise RMS scaling."""
    t0, t1 = _validate_span(t_span, h_initial)
    if rtol <= 0 or h_min <= 0:
        raise ValueError("rtol and h_min must be positive.")
    if h_max is not None and h_max <= 0:
        raise ValueError("h_max must be positive when provided.")

    y_current = np.asarray(y0, dtype=float)
    atol_array = np.asarray(atol, dtype=float) + np.zeros_like(y_current)
    if (atol_array < 0).any():
        raise ValueError("atol must be nonnegative.")

    h = min(h_initial, t1 - t0)
    if h_max is not None:
        h = min(h, h_max)

    t_values = [t0]
    y_values = [y_current.copy()]
    accepted_errors = []
    t_current = t0
    accepted = 0
    rejected = 0
    f_evals = 0
    trials = 0

    while t_current < t1:
        trials += 1
        if trials > max_trials:
            raise RuntimeError("adaptive RK4 exceeded the trial limit.")

        h = min(h, t1 - t_current)
        if h_max is not None:
            h = min(h, h_max)
        if h < h_min:
            raise RuntimeError("adaptive RK4 reached h_min before completing the interval.")

        y_full = rk4_step(f, t_current, y_current, h)
        y_half = rk4_step(f, t_current, y_current, 0.5 * h)
        y_two_half = rk4_step(f, t_current + 0.5 * h, y_half, 0.5 * h)
        f_evals += 12

        scale = atol_array + rtol * np.maximum(np.abs(y_current), np.abs(y_two_half))
        if np.any(scale <= 0):
            raise ValueError("adaptive scale contains non-positive entries.")
        component_error = (y_two_half - y_full) / (15.0 * scale)
        E = float(np.sqrt(np.mean(component_error**2)))

        if E <= 1.0:
            previous_t = t_current
            t_current += h
            if t_current <= previous_t:
                raise RuntimeError("adaptive RK4 failed to advance time.")
            y_current = y_two_half
            t_values.append(t_current)
            y_values.append(y_current.copy())
            accepted_errors.append(E)
            accepted += 1
        else:
            rejected += 1

        factor = max_factor if E == 0.0 else safety * E ** (-0.2)
        factor = min(max_factor, max(min_factor, factor))
        h *= factor

    metadata = {
        "max_estimated_local_error": float(max(accepted_errors) if accepted_errors else 0.0),
        "mean_estimated_local_error": float(np.mean(accepted_errors) if accepted_errors else 0.0),
    }
    return SolverResult(
        np.asarray(t_values),
        np.asarray(y_values),
        "adaptive_rk4",
        SolverStats(accepted, f_evals, rejected_steps=rejected, metadata=metadata),
    )
