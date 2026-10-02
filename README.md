# MTH321 · Week 4 GitHub 文件夹

本包先依据 Tutorial 4 核实任务，再延续 Mathematical Theory（A、B）的项目分工：复现报告数字、校对模型与稳定性公式、整理答辩理解，并为团队最终提交提供可运行证据。

第四周确实有任务，主题为 **Reproduction check, second-wave review, submission**。它不要求凭空引入新模型；重点是用确切代码和版本复现至少一个报告数字到 3 位有效数字，完成提交清单和团队/个人记录。本包还补齐 Project Brief 的匹配精度成本实验。

## 先读这些文件

1. [Tutorial 4 任务与项目结合](notes/00_Week4_Tutorial任务与分工.md)
2. [数学校对、推导与大四学生理解](notes/01_Week4_数学校对与推导.md)
3. [答辩准备与个人复核](notes/02_Week4_答辩与团队复核.md)
4. [实际复现结果与匹配精度成本](notes/04_Week4_实际运行结果.md)
5. [逐项提交检查表](submission/Week4_提交检查表.md)
6. [英文报告结果段落](report/week4_verified_results.md)
7. [独立 LaTeX 数学附录与伪代码](report/week4_math_appendix.tex)

## 如何运行

依赖为 Python 3.10+、NumPy、SciPy 和 Matplotlib。在本文件夹根目录执行：

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python code/run_all.py
```

已装好依赖时：

```text
python code/run_all.py
```

该入口重新生成全部 12 张图、CSV、JSON 和结果页；`code/` 包含前三周需要复用的代码，使第四周包能独立运行。只运行第四周复现与成本实验：

```text
python code/week4_reproduction.py
```

核心输出包括：

```text
REPRODUCED rk4_endpoint_order: report=4.04; run=4.04; units=dimensionless; PASS
REPRODUCED euler_stability_ms: report=0.286; run=0.286; units=ms; PASS
```

这里将 Tutorial 的“3 digits”明确操作为 3 位有效数字；原始数值同时保存在 JSON。图中 4.038 ± 0.009 使用小数位显示，不与 4.04 的复现判据混淆。

## 真实 commit 的复现

本包先记录源码 SHA256 和运行环境；文件打包不会自动产生学生的 Git 贡献记录。真实提交后在干净的 checkout 中运行（把 YOUR_REAL_COMMIT 换成实际 commit）：

```text
python code/week4_reproduction.py --require-clean-commit --expected-commit YOUR_REAL_COMMIT
```

严格模式会拒绝缺失 commit、HEAD 不匹配、源码未跟踪或文件夹有修改的情况。详情见 [团队复核说明](notes/02_Week4_答辩与团队复核.md)。运行证据文件被 `.gitignore` 忽略，以便在干净 commit 上生成；生成后可另存作课程复现记录。

## 每张图的来源

| 图片 | 生成函数/脚本 | 用途 |
|---|---|---|
| `pmsm_currents.png` | `run_all.py:run_project` | 三方法电流与解析参考 |
| `pmsm_convergence.png` | `run_all.py:run_project` | 全网格误差和三方法观测阶 |
| `pmsm_stability_regions.png` | `run_all.py:run_project` | 三种稳定域和 hλ 位置 |
| `pmsm_stability_sweep.png` | `run_all.py:run_project` | 稳定边界两侧的扰动增长/衰减 |
| `pmsm_adaptive.png` | `run_all.py:run_project` | 变化步长、接受/拒绝和局部判据 |
| `synthetic_stiff_case.png` | `run_all.py:run_project` | 人为刚性对照，非实测电机 |
| `lab_logistic.png` | `course_labs.py:run_labs` | Euler 解析特例 |
| `lab_vanderpol.png` | `course_labs.py:run_labs` | 课堂非线性轨迹/失败例 |
| `lab_robertson.png` | `course_labs.py:run_labs` | 参考解与质量守恒 |
| `lab_rc_newton.png` | `course_labs.py:run_labs` | 非线性 Newton/阻尼 |
| `week3_rk4_quantitative.png` | `week3_review.py:run_week3` | Tutorial 3 定量修复及 Week 4 复现对象 |
| `week4_cost_accuracy.png` | `week4_reproduction.py:cost_at_matched_accuracy` | 共同误差上限下的 RHS 成本 |

## 上传与正式提交的边界

可以把本文件夹内容上传到团队仓库；合并同名代码前先保留队友修改。本包不包含课程原始 PDF、虚拟环境或依赖安装目录。

`parameters.json` 仍是教学示例参数；没有使用 `PMSM_Project1_Report.pdf` 的参数。数学附录与英文段落可合入团队报告，它们本身不是已完成的 ≥8 页全队最终报告。成员身份、真实 ICS、最终报告 PDF、展示 PDF、真实 commit 和课程平台记录按检查表完成。由你上传 GitHub。

数学附录已在内置 LaTeX 编辑器打开；编译器返回平台环境错误 `Unable to find standard directories for platform`，因此当前未确认 PDF 编译与排版。源文件保留，可在团队正常 LaTeX 环境中编译并核对。
