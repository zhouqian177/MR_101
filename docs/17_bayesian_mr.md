# 贝叶斯 MR 报告（17_bayesian_mr）

> 对应模块：进阶 3 —— 贝叶斯 MR：cML 约束最大似然法（Xue et al. 2021, AJHG）。
> 脚本：`scripts/32_bayesian_mr.R`；产物：`02.analysis/bayesian/`

## 1. 方法学

**cML（constrained maximum likelihood）** 通过约束最大似然同时识别
"有效工具变量"集合与无效工具变量（存在相关/不相关水平多效性者），
在**多效性普遍存在**时仍可一致估计因果效应：

- 用 BIC/AIC 在候选"无效 IV 数 K"集合间做模型选择
- 模型平均（MA）版本对全部候选模型加权，更稳健
- 无需像 MR-Egger 那样假设"多效性与工具强度无关"

> 实现：MendelianRandomization 包内置 `mr_cML()`（本仓库无需额外安装；
> 独立包 MRcML 亦可从 GitHub xue-hr/MRcML 安装）。

## 2. 数据

LDL-C → CHD（OpenGWAS 在线 harmonised 数据，28 个工具变量）。

## 3. 结果（与 IVW/Egger 对比）

| 方法 | beta | SE | P | OR |
|---|---|---|---|---|
| IVW | 0.403 | 0.220 | 0.067 | 1.496 |
| MR-Egger | 1.040 | 0.426 | 0.015 | 2.831 |
| **cML-BIC(MA)** | **0.467** | **0.149** | **0.0017** | **1.596** |

- cML 未检出无效工具变量（BIC_invalid 为空）
- cML 估计介于 IVW 与 Egger 之间，SE 最小、P 最显著（0.0017）
- 说明：在存在异质性的数据中，cML 通过模型选择给出更精确、稳健的估计

## 4. 结论与教学要点

1. cML 对**相关与不相关多效性均稳健**，优于假设多效性独立的 Egger
2. 当 IVW 与 Egger 结论分歧时，cML 可作为仲裁方法
3. 与 MRMix（扩展 5）互补：MRMix 假设两类 IV 的混合，cML 用约束选择
4. 本结果进一步支持 LDL-C → CHD 的正向因果效应（cML P=0.0017）

## 5. 可复现

```bash
Rscript scripts/32_bayesian_mr.R
```
