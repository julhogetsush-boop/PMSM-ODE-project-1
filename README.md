# MTH321 · PMSM 数值分析项目（Week 1–4）

这个文件夹整合第一到第四周的课程笔记、PMSM 数学推导、Python代码、验证实验、图表与数据、英文报告素材和最终提交检查。全部内容共用一套模型、参数和运行入口。

项目问题是：**在恒速 PMSM 电流初值问题中，方法阶数、稳定域和步长如何共同影响误差与计算成本？怎样用独立参考和可复现代码证明结果可信？**

## 项目模型与参数

研究 d-q 两轴电流状态 $y=(i_d,i_q)^T$。电角速度、输入电压和电机参数在本算例中恒定：

$$\dot i_d=(u_d-Ri_d+\omega_e L_qi_q)/L_d,$$
$$\dot i_q=(u_q-Ri_q-\omega_e L_di_d-\omega_e\psi_f)/L_q.$$

`parameters.json` 是唯一 PMSM 基准参数入口，单位全部为 SI。

| 参数 | 值 | 单位 |
|---|---:|---|
| R | 0.5 | Ω |
| Ld / Lq | 0.003 / 0.004 | H |
| 永磁磁链 ψf | 0.05 | Wb |
| 电角速度 ωe | 1000 | rad/s |
| ud / uq | 0 / 60 | V |
| 初始 id / iq | 0 / 0 | A |
| 末时刻 T | 0.03 | s |

这些是教学示例值。当前模型是外部给定速度的电流子系统；最终项目若使用其他参数，需先修改参数文件，重新运行并同步报告数字。

## 文件夹结构

```text
MTH321_Weeks1-4/
├── README.md                  整个项目的说明与阅读入口
├── AI_Transparency_Log.md     四周统一的AI辅助与人工复核记录
├── parameters.json            PMSM示例参数
├── requirements.txt           Python依赖
├── .gitignore                 环境、缓存和编译中间文件忽略规则
├── code/                      模型、求解器、验证与统一运行入口
├── notes/                     四周笔记、推导、理解与实际结果
├── figures/                   12张由代码生成的PNG
├── results/                   CSV、JSON和实际运行日志
├── report/                    英文结果段落、数学附录与报告数值主张
└── submission/                全项目提交检查表和通用ICS模板
```

## 四周内容与阅读顺序

| 顺序 | 文件 | 内容 |
|---|---|---|
| 0 | [资料依据与四周学习路线](notes/00_资料依据与四周学习路线.md) | 课程来源、各周重点、讲义细节修正 |
| 1 | [Week 1：课程笔记与PMSM建模](notes/01_Week1_课程笔记与PMSM建模.md) | IVP、单位、建模、平衡点、Jacobian、解析解、Euler |
| 2 | [Week 2：数值方法与推导](notes/02_Week2_数值方法与推导.md) | RK4阶条件、隐式Euler、稳定域、复谱步长界、Newton、自适应 |
| 3 | [项目理解与手算示例](notes/03_项目理解与手算示例.md) | 大四学生层面的解释、两步手算、代数与时间误差区别 |
| 4 | [Week 3：Tutorial与定量修复](notes/04_Week3_Tutorial与定量修复.md) | before/after、拟合阶、标准误、理论参考线、浮点误差区 |
| 5 | [Week 4：Tutorial任务与复现](notes/05_Week4_Tutorial任务与复现.md) | 确切版本复现、报告校对、成本证据、课程交付 |
| 6 | [Week 4：数学校对与推导](notes/06_Week4_数学校对与推导.md) | Lyapunov证明、三方法公式核对、误差和结果解释 |
| 7 | [答辩与项目复核](notes/07_答辩与项目复核.md) | 白板问题、真实commit复现、实际复核记录 |
| 8 | [基础实验运行结果](notes/08_基础实验运行结果.md) | PMSM与前两周课堂例题的实际数值和10张图 |
| 9 | [复现与成本运行结果](notes/09_复现与成本运行结果.md) | 报告数字核验和共同误差上限的成本表 |

英文报告素材在 [第三周定量结果](report/week3_quantitative_result.md)、[复现与成本结果](report/week4_verified_results.md) 和 [数学附录LaTeX源文件](report/week4_math_appendix.tex)。完整提交要求见 [统一检查表](submission/提交检查表.md)，个人贡献使用 [通用ICS模板](submission/ICS_贡献填写模板.md)。

## 安装与运行整个项目

需要 Python 3.10+、NumPy、SciPy 和 Matplotlib。在本README所在的根目录打开终端。

Windows PowerShell：

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python code/run_all.py
```

macOS / Linux：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python code/run_all.py
```

已有依赖时直接运行：

```text
python code/run_all.py
```

统一入口依次运行公式/边界检查、PMSM三方法实验、课堂例题、第三周定量修复和第四周复现/成本实验，重写全部12张图、数据和两份实际结果页。成功时最终打印 `PASS. Weeks 1-4 figures, tables, summary.json and result notes regenerated.`

可独立运行的入口：

| 命令 | 用途 |
|---|---|
| `python code/checks.py` | 模型、求解器和独立参考检查 |
| `python code/week3_review.py` | 重做RK4定量收敛图、CSV和JSON |
| `python code/week4_reproduction.py` | 重做报告数字核验、数学恒等式和成本实验 |

程序基于自身文件位置寻找参数和输出；运行不依赖原始课程PDF、MATLAB或Jupyter。所有计算为确定性实验，无随机/Monte Carlo元素。

## 代码文件说明

| 文件 | 负责的内容 |
|---|---|
| `code/model.py` | PMSM的f、Jacobian、平衡点、矩阵指数参考、稳定函数与边界 |
| `code/methods.py` | 通用Euler/RK4/隐式Euler、Newton、步长加倍和差分Jacobian |
| `code/checks.py` | 公式代回、解析特例、非自治阶段时间、拒绝/失败处理、Radau交叉验证 |
| `code/course_labs.py` | Logistic、Van der Pol、Robertson、RC二极管课堂实验 |
| `code/week3_review.py` | RK4细化表、拟合阶与标准误、理论线和舍入区 |
| `code/week4_reproduction.py` | 报告主张比对、Git/源码证据、Lyapunov核对与成本统计 |
| `code/run_all.py` | 整个项目的统一复现入口 |

时间方法和控制器由本项目直接实现。NumPy用于数组与线性求解；SciPy矩阵指数、Radau/DOP853用作参考；Matplotlib负责绘图。`f(t,y)`和`jac(t,y)`通过参数传入通用方法，可以替换模型。

## 图表和数据在哪里生成

| 图文件 | 生成位置 | 对应数据/用途 |
|---|---|---|
| `pmsm_currents.png` | `run_all.py:run_project` | 三方法电流与解析参考 |
| `pmsm_convergence.png` | 同上 | `pmsm_convergence.csv`：全网格误差与三方法阶 |
| `pmsm_stability_regions.png` | 同上 | 三种稳定域与hλ叠图 |
| `pmsm_stability_sweep.png` | 同上 | `pmsm_stability_sweep.csv`：边界两侧扰动 |
| `pmsm_adaptive.png` | 同上 | `pmsm_adaptive_trials.csv`、`pmsm_adaptive_summary.csv` |
| `synthetic_stiff_case.png` | 同上 | `synthetic_stiff_case.csv`：明确标注的人工刚性对照 |
| `lab_logistic.png` | `course_labs.py:run_labs` | `lab_logistic.csv`：解析校准 |
| `lab_vanderpol.png` | 同上 | `lab_vanderpol_week1.csv`及脚本中的失稳实验 |
| `lab_robertson.png` | 同上 | `lab_robertson.csv`：参考解与质量守恒 |
| `lab_rc_newton.png` | 同上 | `lab_rc_trace_full.csv`、`lab_rc_trace_damped.csv`、`lab_rc_step_sweep.csv` |
| `week3_rk4_quantitative.png` | `week3_review.py` | `week3_rk4_convergence.csv`、`week3_review.json` |
| `week4_cost_accuracy.png` | `week4_reproduction.py` | `week4_cost_sweep.csv`、`week4_matched_accuracy.csv` |

`results/summary.json` 汇总四周运行结果与软件版本；`week4_reproduction.json`记录命令、版本、Git状态、源码SHA256和数字核验；`week4_reproduction_line.txt`保存可引用输出行；`week4_run_log.txt`是本次完整运行日志。`pmsm_newton_trace.csv`保留电机首步求根记录，`lab_scalar_stability.csv`保存课堂标量稳定性实验。

## 已验证的主要结果

- 全网格观测阶约为Euler 1.0275、RK4 4.0018、隐式Euler 0.9732。
- 第三周endpoint误差拟合给RK4阶数4.0376，标准误0.0091。
- Euler与RK4稳定性边界分别约为0.285714 ms和2.929622 ms；扰动实验验证指定边界两侧的增长/衰减。
- 自适应控制实际改变步长、记录拒绝，并满足接受步的局部判据。
- 共同误差上限0.005 A下，RK4/Euler/隐式Euler的RHS次数分别为240/30,720/61,440。

全网格误差与endpoint误差是不同指标。局部容差不等同于全局误差上界。成本结论只比较RHS次数，不包含Jacobian/线性求解等全部运行时间。数据和适用范围详见两份结果页。

## 从真实Git版本复现报告数字

`report/report_claims.json`保存报告草稿主张。脚本重新计算后比对3位有效数字，输出：

```text
REPRODUCED rk4_endpoint_order: report=4.04; run=4.04; units=dimensionless; PASS
REPRODUCED euler_stability_ms: report=0.286; run=0.286; units=ms; PASS
```

形成实际提交后，在干净checkout中运行（替换YOUR_REAL_COMMIT）：

```text
python code/week4_reproduction.py --require-clean-commit --expected-commit YOUR_REAL_COMMIT
```

严格模式验证HEAD、源码跟踪与干净状态。尚无真实提交时，普通运行明确记录 `PENDING_REAL_COMMIT`；源码哈希补充记录内容，不替代提交历史。

## GitHub上传和最终提交

将本文件夹内容作为一个仓库根目录，或合入既有仓库的相应目录。README、参数、代码、笔记、图表、CSV/JSON和报告素材一起上传；原始教材、虚拟环境、依赖安装目录和缓存无需加入。生成的复现证据文件被`.gitignore`忽略，可另存作课程记录。

填写真实的 [AI辅助与人工复核记录](AI_Transparency_Log.md) 和ICS，并按统一检查表核对正式报告PDF、展示PDF、真实Git记录与课程提交状态。本项目文件夹提供四周学习和实验的完整材料；正式报告和展示文件仍需据最终项目内容完成。

`report/week4_math_appendix.tex`为可编辑数学附录。内置编译器此前返回平台环境错误，当前未确认PDF编译/排版；请在可用的LaTeX环境中验证。修改参数后，重新审查所有图、报告段落、附录数值表及`report_claims.json`，避免沿用旧数字。
