# MVMR 方法补充报告（29_mvmr_methods）

> 对应模块：拓展 3 —— MR 多变量方法补充（MVMedian/MaxLik）。
> 脚本：`scripts/44_mvmr_methods.R`；产物：`02.analysis/mvmr_supp/`

## 1. 方法学

补充 MendelianRandomization 包中尚未演示的三种方法：

| 方法 | 说明 | 适用场景 |
|---|---|---|
| MV-IVW | 多变量逆方差加权（参考） | 标准多变量 MR |
| MV-Median | 多变量加权中位数 | 对离群工具变量稳健 |
| MaxLik | 最大似然 MR | 效率高，对非正态误差稳健 |

## 2. 结果（模拟数据：X1→Y=0.4, X2→Y=-0.3）

| 方法 | X1 估计 | X2 估计 |
|---|---|---|
| MV-IVW | 0.371 | -0.272 |
| MV-Median | 0.389 | -0.244 |
| MaxLik (X1) | 0.360 | — |

## 3. 解读

- MV-IVW 与 MV-Median 结果一致，相互验证
- MV-Median 对离群工具变量更稳健，适合存在多效性时使用
- 三种方法均接近真实值（0.4, -0.3），验证了方法的有效性

## 4. 可复现

```bash
Rscript scripts/44_mvmr_methods.R
```