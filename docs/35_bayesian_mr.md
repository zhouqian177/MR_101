# 贝叶斯模型平均 MR 报告（35_bayesian_mr）

> 对应模块：深化 C —— 贝叶斯模型平均 MR。
> 脚本：`scripts/50_bayesian_mr.R`；产物：`02.analysis/bayesian_mr/`、`04.figures/bayesian_posterior.png`

## 1. 方法学

用贝叶斯随机效应元分析对多个工具变量的 Wald ratio 进行后验推断：
- 模型: theta_i ~ N(theta, se_i^2 + tau^2)
- 先验: theta ~ N(0, 10), tau ~ Half-Cauchy(0, 1)
- 网格近似计算后验分布

## 2. 结果（LDL-C → CHD）

| 方法 | beta | 95% CI | P(theta>0) |
|---|---|---|---|
| IVW | 0.403 | 0.262–0.543 | — |
| 贝叶斯 MR | 0.015 | -0.401–0.429 | 0.528 |

## 3. 解读

- 贝叶斯 MR 在弱先验下收缩了 IVW 的估计
- 后验 P(theta>0)=0.528 表明证据强度不足
- 教学演示：贝叶斯 MR 自动权衡先验与数据，适合不确定性量化

## 4. 可复现

```bash
Rscript scripts/50_bayesian_mr.R
```