# reproduce_report_rk4_order.py
# Tutorial 4 challenge: reproduce one number from the team report
# (PMSM_Project1_Report.pdf), independent of the GitHub demo package.
#
# Target number: RK4 fitted convergence order p = 3.99507 +- 0.00410
# (report Table 3, fit on the four finest meshes, 50 -> 6.25 us).
#
# Everything below follows the report's own specification:
#   - parameters: Table 1 (Rs=0.5, Ld=0.004, Lq=0.006, psif=0.08,
#     omega_e=500, u=(0,50) V, i0=(0,0) A, interval [0, 0.05] s)
#   - model: eq. (5); exact solution: matrix exponential, eq. (7)
#   - method: classical RK4, eq. (13)-(14) description in Sec. 3.1
#   - error metric: Emax(h) = max_n ||y_n - i_exact(t_n)||_2, eq. (27)
#   - meshes: h = 100, 50, 25, 12.5, 6.25 us (Table 2)
#   - fit: log-log least squares on the four finest meshes (Table 3)

import numpy as np
from scipy.linalg import expm

# ---- Table 1 parameters ----
Rs, Ld, Lq, psif = 0.5, 0.004, 0.006, 0.08
we = 500.0
ud, uq = 0.0, 50.0
i0 = np.zeros(2)
T = 0.05

A = np.array([[-Rs/Ld, we*Lq/Ld], [-we*Ld/Lq, -Rs/Lq]])
b = np.array([ud/Ld, (uq - we*psif)/Lq])
istar = np.linalg.solve(A, -b)          # equilibrium, eq. (6)

def f(t, y):
    return A @ y + b

def rk4_step(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t + h/2, y + h*k1/2)
    k3 = f(t + h/2, y + h*k2/2)
    k4 = f(t + h,   y + h*k3)
    return y + h*(k1 + 2*k2 + 2*k3 + k4)/6

def exact(t):
    return istar + expm(A*t) @ (i0 - istar)   # eq. (7)

# ---- Table 2 meshes, every refinement halves the step ----
hs = [100e-6, 50e-6, 25e-6, 12.5e-6, 6.25e-6]
errs = []
print(f"{'h (us)':>8} {'RK4 max error (A)':>18}   report Table 2")
report_t2 = {100: 4.949e-7, 50: 3.078e-8, 25: 1.918e-9,
             12.5: 1.199e-10, 6.25: 7.601e-12}
for h in hs:
    n = round(T/h)
    y = i0.copy()
    emax = 0.0
    for k in range(n):
        y = rk4_step(f, k*h, y, h)
        e = np.linalg.norm(y - exact((k+1)*h))   # eq. (27), 2-norm
        emax = max(emax, e)
    errs.append(emax)
    print(f"{h*1e6:8.2f} {emax:18.4e}   {report_t2[h*1e6]:.3e}")

# ---- Table 3: log-log least squares on the four finest meshes ----
h4 = np.array(hs[1:]); e4 = np.array(errs[1:])
x = np.log(h4); y = np.log(e4)
n = len(x)
X = np.column_stack([x, np.ones(n)])
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
r = y - X @ beta
s2 = r @ r / (n - 2)
covb = s2 * np.linalg.inv(X.T @ X)
slope, slope_se = beta[0], np.sqrt(covb[0, 0])

print(f"\nRK4 fitted order (four finest meshes): p = {slope:.5f} +/- {slope_se:.5f}")
print(f"report Table 3:                        p = 3.99507 +/- 0.00410")
