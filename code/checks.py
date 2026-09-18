"""Run meaningful checks without a test framework: python code/checks.py."""
import numpy as np
from scipy.integrate import solve_ivp
from methods import (euler_step, rk4_step, implicit_step, solve_fixed,
                     finite_difference_jacobian, adaptive_implicit, NewtonFailure)
from model import load_parameters, matrices, make_model, equilibrium, exact_solution, stability_limit


def run_checks():
    p = load_parameters()
    f, jac = make_model(p)
    A, b = matrices(p)
    assert np.linalg.norm(f(0, equilibrium(p)), np.inf) < 1e-9
    for y in [np.zeros(2), np.array([1.2, -0.7]), equilibrium(p)]:
        diff = finite_difference_jacobian(f, 0, y)
        assert np.linalg.norm(diff-A, np.inf)/max(1, np.linalg.norm(A, np.inf)) < 1e-8
        h = 0.001
        F = lambda w: w-y-h*f(h, w)
        assert np.allclose(finite_difference_jacobian(lambda t, w: F(w), h, y), np.eye(2)-h*A,
                           rtol=1e-8, atol=1e-8)
    # Independent direct solve and analytic first-derivative check.
    y = np.array(p["y0"], dtype=float)
    step = min(0.0001, p["T"]/10)
    got, hist = implicit_step(f, jac, 0, y, step)
    expected = np.linalg.solve(np.eye(2)-step*A, y+step*b)
    assert np.allclose(got, expected, rtol=1e-10, atol=1e-10)
    assert len(hist)-1 <= 2  # normally one; second allowed for roundoff
    z = -0.2
    scalar = lambda t, y: -2*y
    assert np.allclose(euler_step(scalar, 0, np.ones(1), 0.1), [1+z])
    assert np.allclose(rk4_step(scalar, 0, np.ones(1), 0.1), [1+z+z*z/2+z**3/6+z**4/24])
    # Non-autonomous test catches wrong RK stage time / wrong implicit time.
    t, Y, _ = solve_fixed(lambda t, y: np.array([t]), [0], 0, 1, 0.3, "rk4")
    assert t[-1] == 1 and np.all(np.diff(t) > 0)
    assert abs(Y[-1,0]-0.5) < 1e-13
    got, _ = implicit_step(lambda t, y: np.array([t]), lambda t, y: np.zeros((1,1)),
                           0, np.zeros(1), 0.3)
    assert abs(got[0]-0.09) < 1e-13
    # Known decoupled physical special case, independent of matrix exponential.
    zero = dict(p, omega_e=0.0, ud=1.0, uq=2.0)
    t = np.array([0, 0.0001, 0.001, 0.01])
    star = np.array([1.0, 2.0])/zero["R"]
    known = star+(np.zeros(2)-star)*np.exp(-t[:,None]*zero["R"]/np.array([zero["Ld"],zero["Lq"]]))
    assert np.allclose(exact_solution(t, [0,0], zero), known, rtol=1e-12, atol=1e-12)
    # Demonstrate rejection and ensure failed trials do not change initial state.
    t, Y, trials = adaptive_implicit(scalar, lambda t,y: np.array([[-2.0]]), [1],
                                     0, 1, 0.2, atol=0.01, rtol=0)
    assert not trials[0]["accepted"] and trials[1]["t"] == 0
    assert t[-1] == 1 and np.all(np.diff(t) > 0)
    assert all(r["error_ratio"] <= 1 for r in trials if r["accepted"])
    assert abs(Y[-1,0]-np.exp(-2)) < 0.06
    # A failed implicit solve is reported instead of silently accepted.
    try:
        implicit_step(lambda t,y: y+1, lambda t,y: np.ones((1,1)), 0, np.zeros(1), 1)
    except NewtonFailure:
        pass
    else:
        raise AssertionError("Singular implicit equation was silently accepted")
    for bad in [0, -1, np.nan]:
        try:
            solve_fixed(scalar, [1], 0, 1, bad)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid step not rejected")
    # Tight Radau reference versus the exact affine-system formula.
    ts = np.linspace(0, p["T"], 401)
    sol = solve_ivp(f, (0, p["T"]), p["y0"], method="Radau", jac=jac,
                    t_eval=ts, rtol=1e-12, atol=1e-14)
    assert sol.success
    gap = float(np.max(abs(sol.y.T-exact_solution(ts, p["y0"], p))))
    assert gap < 1e-8, gap
    return {"status": "PASS", "radau_vs_matrix_exponential_max_A": gap,
            "checks": "equilibrium; Jf/JF finite differences; direct implicit solve; scalar amplification; stage time; endpoint; decoupled oracle; rejection; failure reporting; invalid h"}


if __name__ == "__main__":
    print(run_checks())
