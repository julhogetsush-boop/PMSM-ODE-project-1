"""PMSM electrical initial value problem.

The model is the two-dimensional linear ODE

    i'(t) = A i(t) + b(t),  i = [i_d, i_q].

For constant voltage and speed, A and b are constant. The exact solution is
therefore available and is used as the oracle for solver validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class PMSMParameters:
    """Motor parameters and baseline voltage input, all in SI units."""

    Rs: float = 0.5
    Ld: float = 4e-3
    Lq: float = 6e-3
    psi_f: float = 0.08
    omega_e: float = 500.0
    pole_pairs: int = 4
    u_d: float = 0.0
    u_q: float = 50.0

    def __post_init__(self) -> None:
        for name in ("Rs", "Ld", "Lq", "psi_f", "omega_e", "u_d", "u_q"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
                raise ValueError(f"{name} must be a finite real scalar.")
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite.")
        if min(self.Rs, self.Ld, self.Lq) <= 0:
            raise ValueError("Rs, Ld and Lq must be strictly positive.")
        if self.psi_f < 0:
            raise ValueError("psi_f must be nonnegative.")
        if (
            isinstance(self.pole_pairs, (bool, np.bool_))
            or not isinstance(self.pole_pairs, Integral)
            or self.pole_pairs <= 0
        ):
            raise ValueError("pole_pairs must be a positive integer.")

    @property
    def omega_m(self) -> float:
        """Mechanical angular speed in rad/s."""
        return self.omega_e / self.pole_pairs


DEFAULT_PARAMETERS = PMSMParameters()
VoltageInput = Union[ArrayLike, Callable[[float], ArrayLike]]


def _pair(values: ArrayLike, name: str) -> NDArray[np.float64]:
    raw = np.asarray(values)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real.")
    pair = np.asarray(values, dtype=float)
    if pair.shape != (2,) or not np.isfinite(pair).all():
        raise ValueError(f"{name} must have shape (2,) and finite entries.")
    return pair.copy()


def _time(t: float, name: str = "t") -> float:
    if not isinstance(t, Real) or not np.isfinite(t):
        raise ValueError(f"{name} must be a finite real scalar.")
    return float(t)


def voltage_input(
    t: float,
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    voltage: Optional[VoltageInput] = None,
) -> NDArray[np.float64]:
    """Return [u_d(t), u_q(t)] in volts."""
    t = _time(t)
    values = (p.u_d, p.u_q) if voltage is None else (
        voltage(t) if callable(voltage) else voltage
    )
    return _pair(values, "voltage")


def jacobian(p: PMSMParameters = DEFAULT_PARAMETERS) -> NDArray[np.float64]:
    """Return A = df/dy for the PMSM current IVP."""
    return np.array(
        [
            [-p.Rs / p.Ld, p.omega_e * p.Lq / p.Ld],
            [-p.omega_e * p.Ld / p.Lq, -p.Rs / p.Lq],
        ],
        dtype=float,
    )


def forcing_vector(
    t: float,
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    voltage: Optional[VoltageInput] = None,
) -> NDArray[np.float64]:
    """Return b(t) in i'(t) = A i(t) + b(t)."""
    u_d, u_q = voltage_input(t, p, voltage=voltage)
    return np.array(
        [
            u_d / p.Ld,
            (u_q - p.omega_e * p.psi_f) / p.Lq,
        ],
        dtype=float,
    )


def rhs(
    t: float,
    y: ArrayLike,
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    voltage: Optional[VoltageInput] = None,
) -> NDArray[np.float64]:
    """Return [di_d/dt, di_q/dt] for y = [i_d, i_q]."""
    i_d, i_q = _pair(y, "y")
    u_d, u_q = voltage_input(t, p, voltage=voltage)
    di_d = (u_d - p.Rs * i_d + p.omega_e * p.Lq * i_q) / p.Ld
    di_q = (u_q - p.Rs * i_q - p.omega_e * (p.Ld * i_d + p.psi_f)) / p.Lq
    return np.array([di_d, di_q], dtype=float)


def equilibrium(
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    voltage: Optional[ArrayLike] = None,
) -> NDArray[np.float64]:
    """Return the constant-input steady current i_star."""
    if callable(voltage):
        raise ValueError("equilibrium requires a constant voltage pair.")
    return np.linalg.solve(jacobian(p), -forcing_vector(0.0, p, voltage=voltage))


def matrix_exponential_2x2(A: ArrayLike, dt: float) -> NDArray[np.float64]:
    """Compute exp(A*dt) for a real 2x2 matrix without scipy."""
    A = np.asarray(A, dtype=float)
    if A.shape != (2, 2):
        raise ValueError("A must have shape (2, 2).")
    dt = _time(dt, "dt")

    identity = np.eye(2)
    mu = 0.5 * np.trace(A)
    B = A - mu * identity
    delta_sq = 0.25 * (A[0, 0] - A[1, 1]) ** 2 + A[0, 1] * A[1, 0]
    delta = np.sqrt(delta_sq + 0j)

    # The closed 2x2 cosh/sinh formula is accurate for the oscillatory PMSM
    # benchmark. In a strong-decay audit, however, cosh and sinh can overflow
    # internally even when the final exp(A t) is small. Use an eigen fallback in
    # that regime so the reference calculation exposes the old failure mode
    # without reproducing it.
    if abs(delta * dt) > 300:
        values, vectors = np.linalg.eig(A)
        exp_A = vectors @ np.diag(np.exp(values * dt)) @ np.linalg.inv(vectors)
    elif abs(delta) < 1e-14:
        exp_A = np.exp(mu * dt) * (identity + dt * B)
    else:
        exp_A = np.exp(mu * dt) * (
            np.cosh(delta * dt) * identity
            + (np.sinh(delta * dt) / delta) * B
        )

    if np.max(np.abs(np.imag(exp_A))) > 1e-10:
        raise FloatingPointError("2x2 matrix exponential produced non-negligible imaginary part.")
    return np.real(exp_A)


def exact_solution(
    t_values: ArrayLike,
    y0: ArrayLike = (0.0, 0.0),
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    t0: float = 0.0,
    voltage: Optional[ArrayLike] = None,
) -> NDArray[np.float64]:
    """Evaluate the constant-input exact formula at scalar or 1-D times."""
    if callable(voltage):
        raise ValueError("exact_solution requires constant voltage and speed.")
    t0 = _time(t0, "t0")
    initial = _pair(y0, "y0")
    raw = np.asarray(t_values)
    if np.iscomplexobj(raw):
        raise ValueError("t_values must be real.")
    times = np.asarray(t_values, dtype=float)
    if times.ndim > 1 or not np.isfinite(times).all():
        raise ValueError("t_values must be a finite scalar or 1-D array.")
    if np.any(times < t0):
        raise ValueError("This IVP interface requires t_values >= t0.")

    A = jacobian(p)
    i_star = equilibrium(p, voltage=voltage)
    offsets = np.atleast_1d(times) - t0
    result = np.empty((offsets.size, 2), dtype=float)
    for k, dt in enumerate(offsets):
        result[k] = i_star + matrix_exponential_2x2(A, float(dt)) @ (initial - i_star)
    return result[0] if times.ndim == 0 else result


def torque(
    y: ArrayLike,
    p: PMSMParameters = DEFAULT_PARAMETERS,
) -> NDArray[np.float64]:
    """Return electromagnetic torque for one or many current states."""
    current = np.asarray(y, dtype=float)
    i_d = current[..., 0]
    i_q = current[..., 1]
    return 1.5 * p.pole_pairs * (p.psi_f * i_q + (p.Ld - p.Lq) * i_d * i_q)


def magnetic_energy(
    y: ArrayLike,
    p: PMSMParameters = DEFAULT_PARAMETERS,
) -> NDArray[np.float64]:
    """Return magnetic co-energy expression used in the power balance."""
    current = np.asarray(y, dtype=float)
    i_d = current[..., 0]
    i_q = current[..., 1]
    return 0.75 * (p.Ld * i_d**2 + p.Lq * i_q**2)


def input_power(
    y: ArrayLike,
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    voltage: Optional[ArrayLike] = None,
) -> NDArray[np.float64]:
    """Return three-phase electrical input power in watts."""
    current = np.asarray(y, dtype=float)
    u_d, u_q = voltage_input(0.0, p, voltage=voltage)
    return 1.5 * (u_d * current[..., 0] + u_q * current[..., 1])


def copper_loss(y: ArrayLike, p: PMSMParameters = DEFAULT_PARAMETERS) -> NDArray[np.float64]:
    """Return stator copper loss in watts."""
    current = np.asarray(y, dtype=float)
    return 1.5 * p.Rs * (current[..., 0] ** 2 + current[..., 1] ** 2)


def flux_error_norm(
    y: ArrayLike,
    reference: ArrayLike,
    p: PMSMParameters = DEFAULT_PARAMETERS,
) -> NDArray[np.float64]:
    """Return ||diag(Ld, Lq) (y-reference)||_2."""
    diff = np.asarray(y, dtype=float) - np.asarray(reference, dtype=float)
    scaled = diff * np.array([p.Ld, p.Lq])
    return np.linalg.norm(scaled, axis=-1)
