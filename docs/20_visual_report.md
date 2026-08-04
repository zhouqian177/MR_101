# MR 结果可视化报告（20_visual_report）

> 对应模块：进阶 6 —— MR 结果可视化报告（批量森林图 + 热图）。
> 脚本：`scripts/35_visual_report.R`；产物：`02.analysis/visual/`、`04.figures/visual_*.png`

## 1. 背景

MR 项目产出大量结果表（多方法 × 多暴露 × 多结局），
需要统一可视化以支持解读与汇报。本模块汇总本仓库各模块结果，
自动生成三类图：

| 图 | 内容 |
|---|---|
| visual_forest_all.png | 各模块 MR 结果森林图（IVW/加权中位数/Egger） |
| visual_heatmap_multi.png | LDL/HDL/TG → CHD 多暴露热图 |
| visual_heatmap_scan.png | LDL-C → 6 结局 MR 扫描热图 |

## 2. 数据来源（自动汇总）

- `02.analysis/mr_results.csv`（两样本 telomere→CHD）
- `02.analysis/opengwas/online/mr_results.csv`（LDL-C→CHD 在线版）
- `02.analysis/mr_scan/scan_results.csv`（MR 扫描）
- `02.analysis/opengwas/multi/comparison_summary.csv`（多暴露对比）

## 3. 实现要点

1. 各模块结果表列结构不同 → **仅保留共同列**（method/b/se/pval/OR 等）后合并
2. 森林图：OR(log 尺度) + 95%CI 误差条 + P<0.05 显著着色
3. 热图：log(OR) 双向配色（蓝=保护、红=危险），标注 OR 与 P 值
4. 结果表导出 `02.analysis/visual/report_summary.csv`

## 4. 产物

- `04.figures/visual_forest_all.png`：22 条结果的统一森林图
- `04.figures/visual_heatmap_multi.png`：LDL OR=1.50 / HDL OR=0.74 / TG OR=1.23 → CHD
- `04.figures/visual_heatmap_scan.png`：LDL-C 对 6 结局的效应矩阵

## 5. 结论与教学要点

1. 可视化是 MR 报告的核心环节：森林图展示效应与精度，热图展示效应谱
2. 热图配色需注意方向语义（log OR 双向色标）
3. 可扩展：批量生成全部暴露×结局组合的热图矩阵（发表级图表）

## 6. 可复现

```bash
Rscript scripts/35_visual_report.R
```
