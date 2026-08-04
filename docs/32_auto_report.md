# MR 结果自动报告报告（32_auto_report）

> 对应模块：深化 3 —— MR 结果自动报告。
> 脚本：`scripts/47_auto_report.R`；产物：`02.analysis/auto_report/`、`04.figures/report_overview.png`

## 1. 功能

自动汇总本仓库所有关键 MR 分析结果，生成论文格式的汇总表与全景图，快速了解仓库全貌。

## 2. 汇总统计

- 仓库分析脚本：34 个
- 仓库文档：29 篇
- 本次汇总分析数：22 条（IVW + 加权中位数）
- 显著结果（IVW）：6 个

## 3. 产出

- `02.analysis/auto_report/paper_ready_table.csv`：论文格式表格（分析/方法/n/beta/SE/OR(95%CI)/P）
- `04.figures/report_overview.png`：仓库结果全景图

## 4. 可复现

```bash
Rscript scripts/47_auto_report.R
```