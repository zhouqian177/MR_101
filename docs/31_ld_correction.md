# MR 的 LD 矩阵校正报告（31_ld_correction）

> 对应模块：深化 2 —— MR 的 LD 矩阵校正。
> 脚本：`scripts/46_ld_correction.R`；产物：`02.analysis/ld_correction/`

## 1. 方法学

标准两样本 MR 假设工具变量间相互独立；当存在 LD 时，IVW 估计的 SE 可能偏小（假阳性）。MendelianRandomization 包的 `mr_input(correlation=)` 参数可接受 LD 矩阵进行校正。

## 2. 模拟结果

| 方法 | beta | SE | P |
|---|---|---|---|
| IVW（未校正 LD） | 0.305 | 0.014 | 2.3e-105 |
| IVW（校正 LD） | 0.307 | 0.016 | 3.0e-81 |

## 3. 解读

- LD 校正后 SE 从 0.014 增大至 0.016（**变化 +15.2%**）
- 校正后 P 值增大，因 LD 使有效工具变量数减少
- 结论：工具间存在 LD 时，不校正会低估 SE，高估显著性

## 4. 可复现

```bash
Rscript scripts/46_ld_correction.R
```