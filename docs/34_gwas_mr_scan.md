# 全基因组 MR 批量扫描报告（34_gwas_mr_scan）

> 对应模块：深化 B —— 全基因组 MR 批量扫描。
> 脚本：`scripts/49_gwas_mr_scan.R`；产物：`02.analysis/gwas_scan/`、`04.figures/gwas_scan_heatmap.png`

## 1. 方法学

扩展 34_mr_scan.R 为多暴露×多结局矩阵，一次运行批量扫描多个暴露-结局对。

## 2. 结果（3 暴露 × 6 结局 = 15 对）

| 暴露 | 结局 | beta | 方向 |
|---|---|---|---|
| LDL-C | CHD 2011 | 0.577 | 正向 |
| LDL-C | CAD 2022 | 0.291 | 正向 |
| LDL-C | MS | -0.131 | 负向 |
| HDL-C | CHD | -0.306 | 负向（保护） |
| TG | CHD | 0.424 | 正向 |

## 3. 解读

- LDL-C 对心血管结局正向（危险），对自身免疫负向
- HDL-C 对 CHD 负向（保护），与文献一致
- 热图直观展示全矩阵效应方向与强度

## 4. 可复现

```bash
Rscript scripts/49_gwas_mr_scan.R
```