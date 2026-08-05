# 药物-靶点网络分析报告（33_drug_target_network）

> 对应模块：深化 A —— 药物-靶点网络分析。
> 脚本：`scripts/48_drug_target_network.R`；产物：`02.analysis/drug_network/`、`04.figures/drug_target_network.png`

## 1. 方法学

整合多个药物靶点蛋白的 pQTL/eQTL 对同一结局的因果效应，构建靶点-疾病效应网络。

## 2. 结果（4 靶点 → CHD）

| 靶点 | 工具数 | beta | 效应方向 | 对应药物 |
|---|---|---|---|---|
| PCSK9 | 1 | 0.013 | 正向（危险） | PCSK9 抑制剂 |
| HMGCR | 2 | -0.040 | 负向（保护） | 他汀类 |
| CETP | 5 | 0.042 | 正向（危险） | CETP 抑制剂 |
| NPC1L1 | 3 | 0.015 | 正向（危险） | 依折麦布 |

## 3. 解读

- HMGCR 表达升高→CHD 风险降低，与他汀抑制 HMGCR 一致
- 网络热图直观展示多靶点×多结局效应矩阵

## 4. 可复现

```bash
Rscript scripts/48_drug_target_network.R
```