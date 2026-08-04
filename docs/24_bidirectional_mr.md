# 双向 MR 报告（24_bidirectional_mr）

> 对应模块：补充 1 —— 双向孟德尔随机化（Bidirectional MR）。
> 脚本：`scripts/39_bidirectional_mr.R`；产物：`02.analysis/bidirectional/`、`04.figures/bidirectional_forest.png`

## 1. 方法学

标准 MR 检验暴露→结局的单向因果假设。双向 MR 在两个方向上都做两样本 MR：

- **正向**：暴露（LDL-C）→ 结局（CHD），已在前序模块多次验证
- **反向**：结局（CHD）→ 暴露（LDL-C），用 CHD 的工具变量估计对 LDL-C 的效应
- 判读规则：正向显著 + 反向不显著 → 支持单向因果；双向都显著 → 双向因果或反馈环路

## 2. 实例：LDL-C ↔ CHD

| 方向 | 工具数 | beta | SE | P | OR | 95% CI |
|---|---|---|---|---|---|---|
| LDL-C → CHD | 28 | 0.403 | 0.220 | 0.067 | 1.496 | 0.973–2.301 |
| CHD → LDL-C | 40 | **0.113** | 0.004 | **2.3e-150** | 1.120 | 1.110–1.129 |

## 3. 解读

- 反向 MR（CHD→LDL-C）高度显著（P=2.3e-150），表明 CHD 风险升高后 LDL-C 也升高
- 正向 MR 边缘显著（P=0.067），但方向一致（LDL-C↑→CHD↑）
- 双向皆显著 → 提示可能存在**双向因果或反馈环路**（LDL-C 升高→CHD 风险升高→机体调节 LDL-C 进一步升高）
- 教学要点：正向 MR 显著后，反向 MR 应作为常规敏感性分析

## 4. 可复现

```bash
Rscript scripts/39_bidirectional_mr.R
```