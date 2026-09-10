"""PMSM electrical initial value problem: parameters, RHS, Jacobian and exact solution.

State and units
---------------
y = [i_d, i_q] in A; time in s; voltage in V; inductance in H.
Amplitude-invariant Park transform, with negative sine in the q row:

    u_d = Rs*i_d + Ld*di_d/dt - omega_e*Lq*i_q
    u_q = Rs*i_q + Lq*di_q/dt + omega_e*(Ld*i_d + psi_f)

Assumptions: balanced three-phase motor, constant parameters and prescribed
constant electrical speed. Saturation, inverter switching and mechanical
acceleration are excluded. Default values are illustrative, not measured data.

Dependencies: numpy, scipy
Install into the same Python interpreter that runs this file:
    python -m pip install -r requirements.txt
Run this file: python model.py

Example: prepare model-independent functions for a numerical integrator::

    p = PMSMParameters()
    f = lambda t, y: rhs(t, y, p)
    jac_f = lambda t, y: jacobian(p)
    reference = exact_solution([0.0, 0.01, 0.05], [0.0, 0.0], p)

The exact solution and equilibrium below require constant voltage. A voltage
function may be supplied to rhs for later numerical experiments. Adding a PI
controller with integral states requires a separate, extended state model.
"""

from dataclasses import dataclass
from numbers import Integral, Real
from typing import Callable, Optional, Union
import os
import site
import sys
import sysconfig


# 某些安装在 Windows 盘符根目录的 Python 将 sys.prefix 记录为 "I:"
# 而不是 "I:\\"，导致 site-packages 被错误地按启动目录解析。
# 仅在该情况下补入此解释器自身的绝对安装路径；保留 -S 的禁用语义。
_drive_root_prefix = os.name == "nt" and len(sys.prefix) == 2 and sys.prefix[1] == ":"
if _drive_root_prefix and not sys.flags.no_site:
    _installed_packages = sysconfig.get_path("purelib")
    if (_installed_packages and os.path.isdir(_installed_packages)
            and _installed_packages not in sys.path):
        site.addsitedir(_installed_packages)

try:
    import numpy as np
    from numpy.typing import ArrayLike, NDArray
    from scipy.linalg import expm
except ModuleNotFoundError as exc:
    # 安装依赖时要使用运行本文件的解释器，避免 pip 安装到另一个 Python。
    # 仅处理缺少顶层依赖；依赖内部的其他导入错误保留原始诊断。
    if exc.name not in {"numpy", "scipy"}:
        raise
    executable = sys.executable.replace("'", "''")
    install_command = f"& '{executable}' -m pip install numpy scipy"
    if _drive_root_prefix:
        # pip 本身也依赖正确的 site-packages 路径，因此临时进入安装根目录。
        install_directory = os.path.dirname(sys.executable).replace("'", "''")
        install_command = (
            f"Push-Location '{install_directory}'\n"
            f"try {{ {install_command} }} finally {{ Pop-Location }}"
        )
    message = (
        f"Missing dependency: {exc.name}\n"
        f"Current Python: {sys.executable}\n\n"
        "Run the following in PowerShell, then run model.py again:\n"
        f"{install_command}"
    )
    if __name__ == "__main__":
        raise SystemExit(message) from None
    raise ModuleNotFoundError(message, name=exc.name) from exc


@dataclass(frozen=True)
class PMSMParameters:
    """Motor parameters and baseline voltage input, all in SI units.

    omega_e is ELECTRICAL angular speed: omega_e = pole_pairs * omega_m.
    u_d and u_q are the constant input used when no voltage override is given.
    """

    Rs: float = 0.5           # 定子每相电阻，ohm
    Ld: float = 4e-3         # d 轴电感，H
    Lq: float = 6e-3         # q 轴电感，H
    psi_f: float = 0.08      # 永磁体磁链，Wb
    omega_e: float = 500.0   # 给定电角速度，rad/s；不要再次乘以极对数
    pole_pairs: int = 4      # 极对数
    u_d: float = 0.0         # 基准 d 轴电压，V
    u_q: float = 50.0        # 基准 q 轴电压，V

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
            raise ValueError("psi_f must be nonnegative in this d-axis convention.")
        if (isinstance(self.pole_pairs, (bool, np.bool_))
                or not isinstance(self.pole_pairs, Integral)
                or self.pole_pairs <= 0):
            raise ValueError("pole_pairs must be a positive integer.")

    @property
    def omega_m(self) -> float:
        """Prescribed mechanical angular speed, in rad/s."""
        return self.omega_e / self.pole_pairs


DEFAULT_PARAMETERS = PMSMParameters()
VoltageInput = Union[ArrayLike, Callable[[float], ArrayLike]]


def _pair(values: ArrayLike, name: str) -> NDArray[np.float64]:
    """Validate one real, finite d-q vector; return an independent array."""
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
    """Return [u_d(t), u_q(t)] in V.

    voltage=None uses p.u_d and p.u_q. An override can be a constant pair or a
    function voltage(t) returning a pair. It must not depend on unknown states.
    """
    t = _time(t)
    values = (p.u_d, p.u_q) if voltage is None else (
        voltage(t) if callable(voltage) else voltage
    )
    return _pair(values, "voltage")


def jacobian(p: PMSMParameters = DEFAULT_PARAMETERS) -> NDArray[np.float64]:
    """Return A = df/dy, shape (2, 2), for prescribed voltage and fixed speed.

    A is independent of time and current for this model. For an integrator
    expecting jac(t, y), use: jac_f = lambda t, y: jacobian(p).
    """
    return np.array([
        [-p.Rs / p.Ld, p.omega_e * p.Lq / p.Ld],
        [-p.omega_e * p.Ld / p.Lq, -p.Rs / p.Lq],
    ], dtype=float)


def forcing_vector(
    t: float,
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    voltage: Optional[VoltageInput] = None,
) -> NDArray[np.float64]:
    """Return b(t) in di/dt = A*i + b(t), in A/s."""
    u_d, u_q = voltage_input(t, p, voltage=voltage)
    return np.array([
        u_d / p.Ld,
        (u_q - p.omega_e * p.psi_f) / p.Lq,
    ])


def rhs(
    t: float,
    y: ArrayLike,
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    voltage: Optional[VoltageInput] = None,
) -> NDArray[np.float64]:
    """Return [di_d/dt, di_q/dt] in A/s for y=[i_d, i_q].

    Suitable for Euler, RK4, implicit Euler or scipy.integrate.solve_ivp.
    Example time-varying drive: voltage=lambda t: [0.0, 50.0 + np.sin(t)].
    """
    i_d, i_q = _pair(y, "y")
    u_d, u_q = voltage_input(t, p, voltage=voltage)

    # 将 d-q 电压方程整理为两个电流导数，保留耦合项和反电动势。
    di_d = (u_d - p.Rs * i_d + p.omega_e * p.Lq * i_q) / p.Ld
    di_q = (u_q - p.Rs * i_q
            - p.omega_e * (p.Ld * i_d + p.psi_f)) / p.Lq
    return np.array([di_d, di_q])


def flux_linkage(
    y: ArrayLike,
    p: PMSMParameters = DEFAULT_PARAMETERS,
) -> NDArray[np.float64]:
    """Return [psi_d, psi_q] in Wb; these are derived quantities, not new states."""
    i_d, i_q = _pair(y, "y")
    return np.array([p.Ld * i_d + p.psi_f, p.Lq * i_q])


def equilibrium(
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    voltage: Optional[ArrayLike] = None,
) -> NDArray[np.float64]:
    """Return constant-input steady current i_star in A, solving A*i_star=-b.

    voltage may override the baseline with a CONSTANT pair. A time-varying
    voltage function has no single equilibrium covered by this function.
    """
    if callable(voltage):
        raise ValueError("equilibrium requires a constant voltage pair.")
    # 求解线性方程组，不显式计算矩阵逆。
    return np.linalg.solve(jacobian(p), -forcing_vector(0.0, p, voltage=voltage))


def exact_solution(
    t_values: ArrayLike,
    y0: ArrayLike = (0.0, 0.0),
    p: PMSMParameters = DEFAULT_PARAMETERS,
    *,
    t0: float = 0.0,
    voltage: Optional[ArrayLike] = None,
) -> NDArray[np.float64]:
    """Evaluate the constant-input exact formula at times t >= t0.

    i(t) = i_star + expm(A*(t-t0)) @ (y0-i_star), with i(t0)=y0.

    Scalar time -> shape (2,). A 1-D time array -> shape (N, 2).
    The expression is analytically exact; expm evaluates it in floating point.
    It works for real or complex-conjugate eigenvalues without branch formulas.
    A time-varying voltage function is rejected; use an independent numerical
    reference or a suitable forcing integral for that case.
    """
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
        raise ValueError("This forward IVP interface requires t_values >= t0.")

    A = jacobian(p)
    i_star = equilibrium(p, voltage=voltage)
    offsets = np.atleast_1d(times) - t0
    result = np.empty((offsets.size, 2), dtype=float)
    for k, dt in enumerate(offsets):
        # 非零初始时刻必须使用 t-t0，而不是直接使用 t。
        result[k] = initial if dt == 0.0 else (
            i_star + expm(A * dt) @ (initial - i_star)
        )
    return result[0] if times.ndim == 0 else result


def _demo() -> None:
    """Print an inspectable baseline example; create no files on import or run."""
    p = PMSMParameters()
    y0 = np.array([0.0, 0.0])
    times = np.array([0.0, 0.005, 0.01, 0.05])
    solution = exact_solution(times, y0, p)
    with np.printoptions(precision=8, suppress=True):
        print("PMSM electrical IVP - illustrative constant-input baseline")
        print("State: [i_d, i_q] (A); time: s")
        print("Parameters:", p)
        print("A = df/dy (1/s):\n", jacobian(p))
        print("b (A/s):", forcing_vector(0.0, p))
        print("Eigenvalues (1/s):", np.linalg.eigvals(jacobian(p)))
        print("Initial derivative (A/s):", rhs(0.0, y0, p))
        print("Equilibrium (A):", equilibrium(p))
        print("Equilibrium residual (A/s):", rhs(0.0, equilibrium(p), p))
        print("\nExact trajectory:")
        print("     t (s)       i_d (A)       i_q (A)")
        for t, y in zip(times, solution):
            print(f"{t:10.6f} {y[0]:13.8f} {y[1]:13.8f}")


if __name__ == "__main__":
    _demo()
