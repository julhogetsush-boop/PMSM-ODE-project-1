# MTH321 · PMSM 项目前两周学习与代码

针对图片中 **Mathematical Theory（A、B）** 的阶段资料：中文课程 notes、PMSM 数学推导、Python 实验、实际运行结果、大四学生层面的理解与答辩准备。

本包依据用户指定的后三项材料：提交清单、学习资料文件夹、分工图片。**没有使用 `PMSM_Project1_Report.pdf` 的参数或结论。`parameters.json` 中是明确标注的教学示例参数，不是实测电机数据。** 模型限定为外部给定恒定电角速度的二状态电流子系统。

## 阅读顺序

1. [资料依据、分工与两周任务](notes/00_阅读依据与两周任务.md)：课程知识如何落到 A、B 的项目交付。
2. [Week 1：课程 notes 与 PMSM 建模](notes/01_Week1_课程笔记与PMSM建模.md)：单位、符号、IVP、平衡点、Jacobian、解析参考解和 Euler。
3. [Week 2：数值方法与数学推导](notes/02_Week2_数值方法与推导.md)：RK4 阶条件、三种稳定域、复特征值步长限制、Newton、步长加倍。
4. [大四学生的理解与答辩准备](notes/03_大四学生理解与答辩.md)：为什么做、结果意味着什么、能推导到什么程度。
5. [实际运行结果](notes/05_运行结果.md)：自动生成的数值、10 张图与解释。
6. [完成情况与 GitHub 上传](notes/04_核对清单与上传说明.md)：你负责上传，包含课堂记录和期末材料的边界。

## 一键运行

需要 Python 3.10+。在本文件所在文件夹打开终端：

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python code/run_all.py
```

如果已在装好依赖的 Python 环境中：

```text
python code/run_all.py
```

仅运行核心检查：`python code/checks.py`。macOS/Linux 可使用 `.venv/bin/python` 替换上面的解释器路径。代码不要求 MATLAB、Jupyter 或原始课程 PDF；没有本机绝对路径。

## 代码与周次

| 文件 | 第一周 | 第二周 |
|---|---|---|
| `code/model.py` | PMSM 的 f、Jacobian、平衡点、参考解 | 稳定函数和步长边界 |
| `code/methods.py` | Euler、统一网格、通用向量接口 | RK4、Newton、隐式 Euler、自适应、差分 Jacobian |
| `code/course_labs.py` | Logistic、标量稳定性、Van der Pol 个体练习 | Van der Pol 失败例、Robertson、RC 二极管阻尼比较 |
| `code/checks.py` | 模型、解析特例与初始数据检查 | Jacobian、阶段时间、拒绝重试、失败报告、独立 Radau 参考 |
| `code/run_all.py` | 统一复现入口 | 生成所有图、数据与结果页 |

`figures/` 保存本次阶段成果的图；`results/` 保存 CSV 和软件版本；运行入口会重写这些已命名输出和 `notes/05_运行结果.md`。不要在自动生成的结果页里保留唯一一份个人批注。

## 参数和结论的范围

`parameters.json` 是唯一 PMSM 基准参数入口。改变参数后重跑所有实验；若收敛阶检查失败，要调整 `run_all.py` 中的步长网格，不能删除失败检查后继续沿用旧图。固定的课堂例题和标为 synthetic 的刚性对照有自己的明确参数。

求解器 Euler/RK4/隐式 Euler/Newton/自适应由本包直接实现；NumPy 用于数组与线性求解，SciPy 的矩阵指数和 Radau/DOP853 用作独立参考，Matplotlib 用于绘图。没有把参考解求解器包装成自己的数值方法。

## 上传前

按 [上传说明](notes/04_核对清单与上传说明.md) 上传**本包内容**，不要把整个原始课程资料、教材或虚拟环境一起上传。填写 [AI 使用与个人贡献记录](AI_Transparency_Log.md) 中真实发生的核验与修改。此包是 AI 辅助整理和实际运行的学习底稿，不自动构成 A、B 已完成的个人贡献证明。
