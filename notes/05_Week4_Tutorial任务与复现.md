# Week 4 · Tutorial 任务与项目复现

## 阅读依据与任务性质

完整阅读 Tutorial 4 的 7 页，并核对 `submission_checklist (1).pdf`、Project Brief 与 Role Description Cards 的相关要求。资料目录只有 Week 1、2 lectures，未发现独立 lecture4；第四周以 Tutorial 4 和期末核对为准。

以下区分课程文件写出的任务与本次整理做出的工作选择：课程要求是依据；用户的请求是查看第四周并按前周形式做成文件夹，由用户上传 GitHub。资料中的论坛发帖、数据库记录、提交和互评是团队/个人应完成的外部任务，不自动代表已执行。

| 阶段 | Tutorial 4 要求 | 在 PMSM 项目中的具体工作 | 本包证据 |
|---|---|---|---|
| 0:00–0:10 | LaTeX 模板、伪代码规范 | 校对状态下标、赋值箭头、一步一行、Newton 容差 | `report/week4_math_appendix.tex` |
| 0:10–0:40 | 挑一个报告数字，从确切 commit 重跑并贴出 3-digit 行 | 选择第三周 RK4 endpoint 拟合阶，同时复核 Euler 稳定界 | `code/week4_reproduction.py`、运行 JSON 与输出行 |
| 0:40–1:00 | 成员比较复现结果，团队同意一个数字；PM 论坛发帖 | 从同一 commit 运行，核对 4.04、误差定义和拟合区间 | 复核表与 PM 草稿，实际发帖待本人执行 |
| 1:00–1:35 | 报告润色、代码注释清理、配对团队交换草稿 | 校对 Model/Methods/Stability、补上成本证据；互评一个优点和一项建议 | 数学附录、英文结果段落、pair-review 记录模板 |
| 1:35–1:50 | 按实时 checklist 逐项检查、口头 wrap-up | 整理最终材料与证据；说明最弱环节 | `submission/提交检查表.md` |

## 本周项目核对内容

核对模型、平衡点、解析参考与 Lyapunov 推导；核对三种稳定域、局部/全局误差、Newton 与自适应估计；从相同代码版本复现报告数字，并能独立解释公式。

| 内容 | 需要核对的事项 | 证据 |
|---|---|---|
| 模型与连续分析 | 电角速度、SI 单位、IVP、A/b、平衡点、矩阵指数、Lyapunov | 数学笔记、公式代回、独立参考 |
| 数值方法 | 三方法阶数、复特征值稳定界、Newton、局部容差、成本指标 | 代码、实验图、CSV/JSON |
| 复现 | 报告数值、源码版本、3位有效数字 | 运行输出、真实commit和源码哈希 |

## 补齐发现的证据缺口

Project Brief 还要求 matched-accuracy cost comparison。前三周的自适应容差表只是容差收紧实验，不能直接当效率排名。本周以 0.01 A 和 0.005 A 为共同误差上限，三种方法沿同一 powers-of-two 网格序列细化，选择第一组达标网格；统计实际 RHS 调用次数，另列 Jacobian 和 Newton 修正次数。

误差定义为各方法自身网格上的最大分量误差，参考为 tight Radau。不同网格的离散最大误差只是采样证据，不声称证明连续时间 supremum。成本不包含 Jacobian/线性求解、参考解和作图；结论限定于 RHS 次数。

## 第四周的课程交付与日期

Tutorial 4 结束要求：至少一个报告数值从确切脚本/命令复现到 3 digits；提交清单完整；PM 发最后 milestone；每人已在 Seminar 2 完成并提交纸质互评和数据库分数/评论。

Project Brief/Checklist 的最终包为：≥8 页英文 LaTeX 报告 PDF（含团队 AI Log 和每人的标注 ICS）、代码 ZIP 加 README、所属展示波次的一套 slides PDF。展示要求 10 分钟加 3 分钟 Q&A；只展示所属波次一次。

资料写的是“最后 tutorial 后的下一个星期一 23:59”，并引用 8 小时宽限及后续每小时扣 10% 的 Syllabus 政策。当前没有最后 tutorial 的实际日期或新版 live form，不能据今天日期推定具体截止日。请以课程实际通知核对。

目录中没有幻灯片引用的 `06_Templates/latex_report_template` 和 live checklist。数学附录是依 Tutorial/伪代码手册编写的独立来源文件，未声称使用缺失模板。
