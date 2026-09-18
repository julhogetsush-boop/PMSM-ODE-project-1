"""恒速 PMSM 电流子系统；SI units; state = [i_d, i_q]."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]


def load_parameters():
    p = json.loads((ROOT / "parameters.json").read_text(encoding="utf-8"))
    values = np.array([p[k] for k in ("R", "Ld", "Lq", "psi_f", "omega_e", "ud", "uq", "T")])
    if not np.all(np.isfinite(values)) or min(p["R"], p["Ld"], p["Lq"], p["T"]) <= 0:
        raise ValueError("Finite parameters and positive R, Ld, Lq, T required")
    if np.shape(p["y0"]) != (2,) or not np.all(np.isfinite(p["y0"])):
        raise ValueError("y0 must contain two finite currents")
    return p


def matrices(p):
    R, Ld, Lq, w = (p[k] for k in ("R", "Ld", "Lq", "omega_e"))
    A = np.array([[-R/Ld, w*Lq/Ld], [-w*Ld/Lq, -R/Lq]])
    b = np.array([p["ud"]/Ld, (p["uq"]-w*p["psi_f"])/Lq])
    return A, b


def make_model(p):
    A, b = matrices(p)
    return (lambda t, y: A @ y + b), (lambda t, y: A.copy())


def equilibrium(p):
    A, b = matrices(p)
    return np.linalg.solve(A, -b)


def exact_solution(t, y0, p):
    """y*= -A^-1 b; y(t)=y*+exp(A*t)(y0-y*). t is elapsed time from 0.

    This analytic formula uses a library matrix exponential, independent of
    our time steppers. No numerical inversion of A is formed.
    """
    A, _ = matrices(p)
    star = equilibrium(p)
    delta = np.asarray(y0, dtype=float) - star
    ts = np.atleast_1d(t)
    Y = np.array([star + expm(A*float(s)) @ delta for s in ts])
    return Y[0] if np.ndim(t) == 0 else Y


def amplification(method, z):
    if method == "euler":
        return 1 + z
    if method == "rk4":
        return 1 + z + z*z/2 + z**3/6 + z**4/24
    if method == "implicit":
        return 1/(1-z)
    raise ValueError(method)


def stability_limit(method, eigenvalues):
    """First positive stability boundary on each stable eigenvalue ray."""
    lam = np.asarray(eigenvalues, dtype=complex)
    if np.any(lam.real >= 0):
        raise ValueError("This routine requires strictly decaying eigenvalues")
    if method == "implicit":
        return np.inf
    if method == "euler":
        return float(np.min(-2*lam.real / abs(lam)**2))
    if method != "rk4":
        raise ValueError(method)
    # Along a left-half-plane ray classical RK4 has the connected interval
    # from h=0 to its first exit. Bracket that exit, then bisect.
    lo, hi = 0.0, 1/np.max(abs(lam))
    while np.max(abs(amplification(method, hi*lam))) <= 1:
        hi *= 2
    for _ in range(70):
        mid = (lo+hi)/2
        if np.max(abs(amplification(method, mid*lam))) <= 1:
            lo = mid
        else:
            hi = mid
    return float((lo+hi)/2)
