"""Regenerate every figure and result: python code/run_all.py.

Paths are relative to this repository, not the current working directory.
All parameters are SI; no random experiments are used.
"""
import json
import platform
from pathlib import Path
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import (load_parameters, matrices, make_model, equilibrium,
                   exact_solution, amplification, stability_limit)
from methods import solve_fixed, adaptive_implicit, implicit_step
from checks import run_checks
from course_labs import run_labs, save_csv, save_figure

ROOT = Path(__file__).resolve().parents[1]
METHODS = ["euler", "rk4", "implicit"]
NAMES = {"euler": "Explicit Euler", "rk4": "RK4", "implicit": "Implicit Euler"}
COLORS = {"euler": "#0072B2", "rk4": "#009E73", "implicit": "#D55E00"}


def run_project(figures, results):
    p = load_parameters()
    f, jac = make_model(p)
    A, b = matrices(p)
    lam = np.linalg.eigvals(A)
    star, y0, T = equilibrium(p), np.array(p["y0"]), p["T"]
    limits = {m: stability_limit(m, lam) for m in METHODS}
    summary = {"parameters": p, "equilibrium_A": star.tolist(),
               "jacobian_per_s": A.tolist(),
               "eigenvalues_per_s": [[float(z.real),float(z.imag)] for z in lam],
               "euler_h_boundary_s": limits["euler"], "rk4_h_boundary_s": limits["rk4"],
               "decay_rate_ratio": float(max(abs(lam.real))/min(abs(lam.real)))}
    # Whole-grid norm avoids hiding transient error at a near-equilibrium endpoint.
    n_values = [240, 480, 960, 1920, 3840]
    conv = []
    for n in n_values:
        h = T/n
        oracle = exact_solution(np.linspace(0,T,n+1), y0, p)
        for m in METHODS:
            t, Y, counts = solve_fixed(f, y0, 0, T, h, m, jac)
            err = abs(Y-oracle)
            conv.append({"method": m, "N": n, "h_s": h,
                         "max_grid_Linf_error_A": float(np.max(err)),
                         "endpoint_Linf_error_A": float(np.max(err[-1])), **counts})
    save_csv(results/"pmsm_convergence.csv", conv)
    fig, ax = plt.subplots(figsize=(7,4.5))
    orders = {}
    for m in METHODS:
        rows = [r for r in conv if r["method"] == m]
        hs = np.array([r["h_s"] for r in rows])
        errs = np.array([r["max_grid_Linf_error_A"] for r in rows])
        order = float(np.polyfit(np.log(hs[-3:]), np.log(errs[-3:]), 1)[0])
        orders[m] = order
        ax.loglog(hs, errs, "o-", color=COLORS[m], label=f"{NAMES[m]}: p={order:.3f}")
    ax.set(xlabel="Uniform step h (s)", ylabel="Max grid/component error (A)",
           title="PMSM: convergence against matrix-exponential reference")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend()
    save_figure(fig, figures/"pmsm_convergence.png")
    summary["observed_orders_last3"] = orders
    assert 0.9 < orders["euler"] < 1.2 and 0.8 < orders["implicit"] < 1.1
    assert 3.8 < orders["rk4"] < 4.2
    # Two-component trajectories.
    fig, axes = plt.subplots(1,2,figsize=(11,4))
    dense = np.linspace(0,T,1201)
    exact = exact_solution(dense,y0,p)
    for j, ax in enumerate(axes):
        ax.plot(dense*1000, exact[:,j], "k--", lw=1.5, label="Matrix exponential")
        for m in METHODS:
            t,Y,_ = solve_fixed(f,y0,0,T,T/600,m,jac)
            ax.plot(t*1000,Y[:,j],color=COLORS[m],label=NAMES[m],lw=1)
        ax.set(xlabel="Time (ms)",ylabel=f"i_{'dq'[j]} (A)",title=f"{'dq'[j]}-axis current; h={T/600*1e6:.1f} us")
        ax.legend(fontsize=8)
    save_figure(fig,figures/"pmsm_currents.png")
    # Stability regions with the same actual eigenvalues and step sizes.
    X,Y = np.meshgrid(np.linspace(-4,2,550),np.linspace(-4,4,550))
    Z = X+1j*Y
    fig,axes = plt.subplots(1,3,figsize=(12,4.4))
    probes = [(0.8*limits["euler"],"o","0.8 h_E"),
              (1.2*limits["euler"],"x","1.2 h_E"),
              (1.2*limits["rk4"],"+","1.2 h_RK4")]
    for ax,m in zip(axes,METHODS):
        with np.errstate(divide="ignore",invalid="ignore"):
            mask = abs(amplification(m,Z)) <= 1
        ax.contourf(X,Y,mask.astype(float),levels=[0.5,1.5],colors=[COLORS[m]],alpha=0.22)
        ax.contour(X,Y,mask.astype(float),levels=[0.5],colors=[COLORS[m]],linewidths=1)
        for h,mark,label in probes:
            ax.scatter((h*lam).real,(h*lam).imag,marker=mark,s=40,label=label)
        ax.axhline(0,color="grey",lw=0.5); ax.axvline(0,color="grey",lw=0.5)
        ax.set(xlabel="Re(h lambda)",ylabel="Im(h lambda)",title=NAMES[m],xlim=(-4,2),ylim=(-4,4))
        ax.set_aspect("equal")
    axes[0].legend(fontsize=8)
    save_figure(fig,figures/"pmsm_stability_regions.png")
    # Actual perturbation experiment across each explicit stability boundary.
    fig,axes = plt.subplots(1,2,figsize=(11,4))
    sweep = []
    for ax,m in zip(axes,["euler","rk4"]):
        for ratio in [0.8,1.2]:
            h = ratio*limits[m]
            horizon = 120*h
            t,Y,_ = solve_fixed(f,star+np.array([0.01,0]),0,horizon,h,m,jac)
            # Weighted flux perturbation norm is compatible with the Lyapunov proof.
            weights = np.array([p["Ld"],p["Lq"]])
            norm = np.linalg.norm((Y-star)*weights,axis=1)
            amp = float(np.max(abs(amplification(m,h*lam))))
            ax.semilogy(t*1000,np.maximum(norm/norm[0],1e-16),label=f"h/hcrit={ratio}; rho={amp:.3f}")
            sweep.append({"method":m,"h_ratio":ratio,"h_s":h,"spectral_radius":amp,
                          "T_s":horizon,"final_relative_flux_perturbation":float(norm[-1]/norm[0])})
        ax.set(xlabel="Time (ms)",ylabel="Weighted perturbation / initial",title=f"{NAMES[m]}: 120 steps")
        ax.legend(fontsize=8)
    save_csv(results/"pmsm_stability_sweep.csv",sweep)
    save_figure(fig,figures/"pmsm_stability_sweep.png")
    for row in sweep:
        if row["h_ratio"] < 1:
            assert row["spectral_radius"] < 1 and row["final_relative_flux_perturbation"] < 1
        else:
            assert row["spectral_radius"] > 1 and row["final_relative_flux_perturbation"] > 1
    # Adaptive implicit Euler; count RHS calls including rejected trials.
    fig,axes = plt.subplots(3,1,figsize=(8,8),sharex=True)
    adaptive_rows = []
    for rtol in [1e-3,1e-4,1e-5]:
        calls = {"rhs":0}
        def counted_rhs(t,y):
            calls["rhs"] += 1
            return f(t,y)
        t,Y,trials = adaptive_implicit(counted_rhs,jac,y0,0,T,T/10,atol=rtol/100,rtol=rtol)
        accepted = [r for r in trials if r["accepted"]]
        exact = exact_solution(t,y0,p)
        row = {"rtol":rtol,"atol_A":rtol/100,"accepted":len(accepted),
               "rejected":len(trials)-len(accepted),"rhs_calls_including_rejects":calls["rhs"],
               "h_min_s":min(r["h"] for r in accepted),"h_max_s":max(r["h"] for r in accepted),
               "max_accepted_local_ratio":max(r["error_ratio"] for r in accepted),
               "max_grid_Linf_global_error_A":float(np.max(abs(Y-exact))),
               "endpoint_Linf_error_A":float(np.max(abs(Y[-1]-exact[-1])))}
        adaptive_rows.append(row)
        if rtol == 1e-4:
            save_csv(results/"pmsm_adaptive_trials.csv",trials)
            axes[0].plot(t*1000,Y[:,0],label="adaptive i_d")
            axes[0].plot(t*1000,exact[:,0],"k--",label="reference")
            axes[0].set(ylabel="Current (A)",title="Adaptive implicit Euler; rtol=1e-4, atol=1e-6 A")
            axes[0].legend()
            axes[1].semilogy(t[1:]*1000,np.diff(t)*1000)
            axes[1].set(ylabel="Accepted h (ms)")
            axes[2].plot([r["t"]*1000 for r in accepted],[r["error_ratio"] for r in accepted],".",label="accepted")
            rejected = [r for r in trials if not r["accepted"] and np.isfinite(r["error_ratio"])]
            axes[2].plot([r["t"]*1000 for r in rejected],[r["error_ratio"] for r in rejected],"x",label="rejected")
            axes[2].axhline(1,color="black",ls="--")
            axes[2].set(yscale="log",ylabel="Scaled local estimate E",xlabel="Time (ms)")
            axes[2].legend()
        assert row["max_accepted_local_ratio"] <= 1 and row["rejected"] > 0
        assert row["h_max_s"] > 1.5*row["h_min_s"]
    save_csv(results/"pmsm_adaptive_summary.csv",adaptive_rows)
    save_figure(fig,figures/"pmsm_adaptive.png")
    summary["adaptive"] = adaptive_rows
    assert adaptive_rows[-1]["max_grid_Linf_global_error_A"] < adaptive_rows[0]["max_grid_Linf_global_error_A"]
    # A clearly labelled synthetic stiff parameter case; not this baseline motor.
    stiff = dict(p,R=0.5,Ld=1e-4,Lq=0.02,omega_e=0.0,ud=1.0,uq=1.0)
    sf,sj = make_model(stiff)
    fig,axes = plt.subplots(1,2,figsize=(11,4))
    stiff_rows = []
    for m,h in [("euler",0.0002),("euler",0.0006),("implicit",0.005)]:
        t,Y,_ = solve_fixed(sf,[0,0],0,0.05,h,m,sj)
        for j,ax in enumerate(axes):
            ax.plot(t*1000,Y[:,j],label=f"{NAMES[m]}, h={h*1000:g} ms")
        stiff_rows.append({"method":m,"h_s":h,"endpoint_error_A":float(np.max(abs(Y[-1]-exact_solution(0.05,[0,0],stiff))))})
    tt = np.linspace(0,0.05,501)
    ref = exact_solution(tt,[0,0],stiff)
    for j,ax in enumerate(axes):
        ax.plot(tt*1000,ref[:,j],"k--",label="Reference")
        ax.set(xlabel="Time (ms)",ylabel=f"i_{'dq'[j]} (A)",title=["Fast mode: lambda=-5000/s","Slow mode: lambda=-25/s"][j])
        ax.legend(fontsize=7)
    axes[0].set_yscale("symlog",linthresh=3)
    save_csv(results/"synthetic_stiff_case.csv",stiff_rows)
    save_figure(fig,figures/"synthetic_stiff_case.png")
    summary["synthetic_stiff_case"] = {"decay_ratio":200,"euler_h_boundary_s":0.0004,
                                        "description":"Deliberate numerical stress test; not measured motor data"}
    _,hist = implicit_step(f,jac,0,y0,0.001)
    save_csv(results/"pmsm_newton_trace.csv",hist)
    summary["pmsm_first_step_newton_updates"] = len(hist)-1
    return summary


def write_result_note(summary):
    p, labs = summary["project"], summary["course_labs"]
    orders = p["observed_orders_last3"]
    lines = ["# 前两周实际运行结果", "", "本页由 `python code/run_all.py` 自动生成。参数见 `parameters.json`，是教学示例参数。",
             "", "## PMSM 基准算例", "",
             f"- 平衡电流 [i_d*, i_q*] = {p['equilibrium_A']} A。",
             f"- 特征值 [实部, 虚部] = {p['eigenvalues_per_s']}，单位 1/s。",
             f"- Euler 稳定性边界：{p['euler_h_boundary_s']*1000:.9f} ms；RK4：{p['rk4_h_boundary_s']*1000:.9f} ms。边界等号不保证衰减。",
             f"- 衰减率之比：{p['decay_rate_ratio']:.4g}；当前基准不具备明显的快慢衰减尺度分离。",
             f"- 网格最大误差的最后三点拟合阶：Euler {orders['euler']:.5f}；RK4 {orders['rk4']:.5f}；隐式 Euler {orders['implicit']:.5f}。",
             f"- 独立 Radau 与矩阵指数参考解的最大差：{summary['checks']['radau_vs_matrix_exponential_max_A']:.3e} A。",
             f"- PMSM 首步 Newton 修正次数：{p['pmsm_first_step_newton_updates']}。这验证线性残差的性质，不能据此声称验证了非线性二次收敛。",
             "", "![电流](../figures/pmsm_currents.png)", "", "两条电流均与矩阵指数解比较；同一步长下可看出一阶方法的相位和阻尼误差。",
             "", "![收敛](../figures/pmsm_convergence.png)", "", "误差取全部网格点、全部分量的最大绝对误差，避免末时刻接近平衡点而掩盖暂态误差。拟合只取最后三点，完整数据保留在 results。",
             "", "![稳定域](../figures/pmsm_stability_regions.png)", "", "阴影满足 |R(z)|≤1；同一特征值在不同步长下的位置说明方法允许的步长不同。隐式 Euler 的完整稳定域为圆盘 |1-z|<1 的外部及边界。",
             "", "![步长扫描](../figures/pmsm_stability_sweep.png)", "", "每条曲线运行 120 步，时长随步长变化；用途是检验扰动的增长或衰减，不是比较相同末时刻的精度。",
             "", "## 自适应隐式 Euler", "", "| rtol | 接受/拒绝 | RHS 次数 | h 最小/最大 (ms) | 最大局部 E | 全程网格误差 (A) |", "|---|---:|---:|---:|---:|---:|"]
    for r in p["adaptive"]:
        lines.append(f"| {r['rtol']:.0e} | {r['accepted']}/{r['rejected']} | {r['rhs_calls_including_rejects']} | {r['h_min_s']*1e3:.5g}/{r['h_max_s']*1e3:.5g} | {r['max_accepted_local_ratio']:.4f} | {r['max_grid_Linf_global_error_A']:.5g} |")
    lines += ["", "每次试算含一次整步和两次半步；RHS 次数包括拒绝步和残差求值，不计 Jacobian、线性分解、参考解或绘图成本。此表展示容差收紧的效果，不是匹配全局精度的效率排名。",
              "", "![自适应](../figures/pmsm_adaptive.png)", "", "接受步的 E≤1，拒绝步保留原时间和原状态后重试。E 是局部误差估计，不能当作全局误差上界。",
              "", "![刚性对照](../figures/synthetic_stiff_case.png)", "", "人为对照参数给出 -5000 与 -25 的两个实特征值，衰减率之比为 200；这与基准的振荡限制不同。隐式大步稳定，但快暂态仍可能不准确。",
              "", "## 课堂例题复现", "",
              f"- Logistic Euler 最后三点阶：{labs['logistic_order_last3']:.5f}；正确手算两步 0.595 和 0.7069195 已通过（讲义第二步有算术笔误，见资料核对页）。",
              f"- Tutorial 1 Van der Pol：mu=1、y0=(0.5,0)、T=10；参考 y1(T)={labs['vdp_week1_y1_T']:.10f}；规定四步长的端点 L2 拟合阶={labs['vdp_week1_order_prescribed4']:.5f}。额外细网格保留在 CSV。",
              f"- Week 2 Van der Pol：mu=10、y0=(2,0)、h=0.1，在 t={labs['vdp_failure_last_t']:.3f} 时触发状态幅值 1e6 的停止阈值；不是参考解发散。",
              f"- Robertson：隐式 h=0.1 到 T=40；最大质量缺陷 {labs['robertson_mass_error']:.3e}，末点分量最大误差 {labs['robertson_endpoint_error_max']:.3e}。质量守恒与轨迹精度单独衡量。",
              f"- RC 二极管首步 h=1 ms：完整 Newton {labs['rc_first_step_full_updates']} 次，阻尼 Newton {labs['rc_first_step_damped_updates']} 次；两者到达同一个隐式根。",
              f"- 在给定离散扫描网格上，每步至多 10 次 Newton 修正的最大步长：完整 Newton {labs['rc_full_hmax_ms']:g} ms，阻尼 Newton {labs['rc_damped_hmax_ms']:g} ms。这不是连续 h 的精确临界值。",
              "", "![Logistic](../figures/lab_logistic.png)", "", "解析解校准 Euler；一阶收敛是可重复的数值证据。",
              "", "![Van der Pol](../figures/lab_vanderpol.png)", "", "左图是第一周指定初值的相图；右图是第二周不同参数的失败例。两者不可混用。",
              "", "![Robertson](../figures/lab_robertson.png)", "", "轨迹比较与质量缺陷并列，说明守恒检查不能代替独立参考解。",
              "", "![RC Newton](../figures/lab_rc_newton.png)", "", "阻尼减少首步过冲；横轴是 Newton 修正次数，不是物理时间。完整迭代表见 results/lab_rc_step_sweep.csv。",
              "", "## 验证范围", "", "所有输出由本次代码实际计算；完整设置和软件版本见 results/summary.json。数值验证只覆盖这些算例，不声称证明任意非线性系统的稳定性或任意参数下的正确性。学生本人应复跑并填写自己的理解、修改和真实贡献记录。"]
    (ROOT/"notes/05_运行结果.md").write_text("\n".join(lines)+"\n",encoding="utf-8")


def main():
    figures,results = ROOT/"figures",ROOT/"results"
    for directory in [figures,results,ROOT/"notes"]:
        directory.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    print("Checking formulas, edge cases and independent oracle...",flush=True)
    checks = run_checks()
    print("Running PMSM experiments...",flush=True)
    project = run_project(figures,results)
    print("Running Week 1-2 classroom examples...",flush=True)
    labs = run_labs(figures,results)
    summary = {"checks":checks,"project":project,"course_labs":labs,
               "versions":{"python":platform.python_version(),"numpy":np.__version__,
                           "scipy":scipy.__version__,"matplotlib":matplotlib.__version__}}
    (results/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False),encoding="utf-8")
    write_result_note(summary)
    print("PASS. All figures, CSV tables, summary.json and the result note regenerated.")
    print("Observed orders:",project["observed_orders_last3"])
    print("Euler / RK4 stability boundaries (ms):",project["euler_h_boundary_s"]*1000,project["rk4_h_boundary_s"]*1000)
    print("RC first-step Newton updates:",labs["rc_first_step_full_updates"],labs["rc_first_step_damped_updates"])


if __name__ == "__main__":
    main()
