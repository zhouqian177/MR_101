# 数据准备报告（02_data_prep）

> 对应模块：00.data + 01_data_prep —— 暴露/结局 GWAS 数据、工具变量提取、LD clumping、harmonise。
> 脚本：`scripts/01_data_prep.R`；产物：`02.analysis/instruments.txt`、`harmonised.csv`、`data_prep_summary.txt`

## 1. 数据来源

| 角色 | 数据集 | 说明 |
|---|---|---|
| 暴露 | telomere_length.txt | 端粒长度 GWAS（Codd et al. 2013，TwoSampleMR 官方示例） |
| 结局 | cardiogram.txt | 冠心病 CHD（CARDIoGRAM 联盟，TwoSampleMR 官方示例） |
| LD 参考 | CEU HapMap3 | 借用 01.GWAS 项目 01.qc/ceu_raw 面板（165 人，145.8 万 SNP） |

> 说明：由于当前环境 OpenGWAS API（api.opengwas.io）被代理阻断，无法在线提取
> 完整 GWAS 汇总数据，故采用 TwoSampleMR 包内置的真实 GWAS 示例数据。
> 在线流程（`extract_instruments` / `extract_outcome_data`）与本地流程完全一致。

## 2. 工具变量提取与 LD clumping

- 候选工具变量：31 个（官方示例已筛选的端粒长度相关 SNP 列表）
- LD clumping 参数：CEU HapMap3 参考面板，r² < 0.001，窗口 10 Mb
- clumping 后保留：**17 个**（去除 LD 冗余 SNP）

## 3. Harmonise

- 统一暴露/结局的效应等位基因与效应方向
- harmonise 后保留：**17 个 SNP**（全部 mr_keep=TRUE）

## 4. 弱工具变量检验（F 统计量）

F = (beta/se)²，规则：F > 10 视为无弱工具变量问题。

| 指标 | 值 |
|---|---|
| 最小 F | 0 |
| 平均 F | 8.6 |
| 最大 F | 29.6 |
| F < 10 的 SNP 数 | 11 / 17 |

**结论**：本示例数据样本量较小（~2 万），平均 F=8.6 偏弱，存在弱工具变量
风险。实际研究建议先过滤 F > 10；本项目为教学目的保留全部 17 个 SNP，
并在分析报告中讨论该局限性（弱工具变量会导致 IVW 估计偏向 0，即偏保守）。

## 5. 产物清单

| 文件 | 内容 |
|---|---|
| 02.analysis/instruments.txt | 17 个工具变量的暴露/结局效应、P 值、F 统计量 |
| 02.analysis/harmonised.csv | harmonise 后完整数据（供后续分析使用） |
| 02.analysis/data_prep_summary.txt | 数据准备摘要 |
| 02.analysis/clump.clumped | LD clumping 结果 |

## 6. 可复现

```bash
Rscript scripts/01_data_prep.R
```
