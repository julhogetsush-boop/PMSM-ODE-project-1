# Week 3 · GitHub 内容清单

## 本周新增并应提交的文件

| 文件 | 用途 |
|---|---|
| `code/week3_review.py` | 从模型重新计算 RK4 细化表、拟合阶与标准误，并生成定量图 |
| `figures/week3_rk4_quantitative.png` | Tutorial 3 的 after 图；关闭 convergence plot feedback |
| `results/week3_rk4_convergence.csv` | 全部原始步长、误差、拟合选择和舍入区标志 |
| `results/week3_review.json` | 图中核心数字与解释范围 |
| `notes/06_Week3_Tutorial与定量修复.md` | Tutorial 要求、before/after、finding/fix plan 和限制 |
| `report/week3_quantitative_result.md` | 可放入英文报告的定量段落与图注草稿 |
| `README.md` | 第三周阅读顺序、文件结构和独立运行方法 |
| `code/model.py`、`code/methods.py`、`parameters.json` | 第三周脚本独立运行所需的既有模型、求解器和示例参数 |

`figures/pmsm_convergence.png` 是 before 图，来自前两周包，保留以便说明修改前后的差别。

## 上传前本地检查

在仓库根目录执行：

```text
python code/week3_review.py
```

成功标准：命令正常结束，打印的 JSON 中 slope 接近 4、slope_standard_error 小于 0.05。随后核对：

- `figures/week3_rk4_quantitative.png` 能打开，图例和 annotation 没有遮挡数据；
- `results/week3_review.json` 的 slope 接近 4，且和报告段落一致；
- `git diff` 中没有 `.venv/`、`__pycache__/`、原始教材或个人敏感信息；
- 提交说明描述真实修改，例如 `Add Tutorial 3 quantitative RK4 review figure`；
- A、B 各自的实际修改和核验由真实 commit/issue/ICS 记录。

## 建议 Git 提交边界

若团队希望一条可审查提交，可将上表文件作为同一 commit：代码、生成数据、图和解释一起出现，审查者能从 finding 追到脚本和输出。不要只上传 PNG，否则 Tutorial 4 无法从确切脚本复现数字。

本包没有推送远程仓库。可直接上传整个 `MTH321_Week3_GitHub` 文件夹，或把第三周新增文件合入团队现有结构；合入前先处理同名文件，保留队友改动。
