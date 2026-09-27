# MTH321 · Week 3 GitHub 文件夹

本文件夹只包含第三周 **Tutorial 3: Scientific visualization and first-wave review** 对应的 GitHub 内容。它把一条具体反馈转成可复现的定量结果：原收敛图缺少理论斜率、拟合不确定性与浮点误差区，因此重新生成一张只验证 RK4 四阶收敛的 after 图。

## 建议阅读顺序

1. `notes/Week3_Tutorial与定量修复.md`：Tutorial 3 要求、before/after、数学解释和 finding/fix plan。
2. `figures/pmsm_convergence.png`：修改前的 before 图。
3. `figures/week3_rk4_quantitative.png`：修改后的定量图。
4. `results/week3_rk4_convergence.csv`：全部步长和误差数据。
5. `results/week3_review.json`：拟合范围、斜率、标准误和舍入区边界。
6. `report/week3_quantitative_result.md`：可合入英文项目报告的结果段落和图注。
7. `notes/Week3_GitHub清单.md`：上传和个人记录检查。

## 文件结构

```text
MTH321_Week3_GitHub/
├── .gitignore
├── README.md
├── AI_Transparency_Log.md
├── parameters.json
├── requirements.txt
├── code/
│   ├── model.py
│   ├── methods.py
│   └── week3_review.py
├── figures/
│   ├── pmsm_convergence.png
│   └── week3_rk4_quantitative.png
├── results/
│   ├── week3_rk4_convergence.csv
│   └── week3_review.json
├── notes/
│   ├── Week3_Tutorial与定量修复.md
│   └── Week3_GitHub清单.md
└── report/
    └── week3_quantitative_result.md
```

`model.py` 和 `methods.py` 是第三周脚本运行所需的最小依赖，并非本周新提出的方法。课程 Tutorial PDF 没有复制到文件夹中。

## 运行第三周代码

需要 Python 3.10+。在本文件夹中执行：

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python code/week3_review.py
```

已有依赖时可直接运行：

```text
python code/week3_review.py
```

脚本会重写第三周 PNG、CSV 和 JSON，并检查拟合阶接近 4、标准误足够小、细网格误差出现非单调浮点平台。当前实际结果为 **4.0376 ± 0.0091**。

## 使用范围

`parameters.json` 使用明确标注的教学示例参数，不是团队最终实测电机参数。真实 Seminar 1 反馈、课程数据库记录、peer-review 表和 Git 历史没有被虚构。上传前请由 A、B 填写本人真正完成的复核与修改。
