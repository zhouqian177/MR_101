# MR 扫描报告（19_mr_scan）

> 对应模块：进阶 5 —— 全基因组 MR 扫描（MR Scan / PheWAS-MR 批量）。
> 脚本：`scripts/34_mr_scan.R`；产物：`02.analysis/mr_scan/`、`04.figures/mr_scan_forest.png`

## 1. 方法学

MR 扫描（MR Scan）用一个暴露的工具变量，对**多个结局数据集**批量执行
两样本 MR，快速评估该暴露的因果效应谱（phenome-wide MR），
常用于：

- 暴露安全性筛查（如药物靶点蛋白对多种疾病的风险）
- 发现新的因果关联（如 LDL-C 与哪些疾病相关）
- 检验特异性（危险因素应只影响相关结局）

## 2. 数据（OpenGWAS 在线）

- 暴露：LDL-C（ieu-b-110），28 个工具变量（LD clumping 后）
- 结局（6 个数据集）：CHD×2（ieu-a-7/ieu-a-8）、CAD（ebi-a-GCST005194）、
  多发性硬化（ieu-b-18）、淋巴白血病（ieu-b-4956）、HOMA-IR（ieu-b-118）、
  雌二醇（ieu-b-4872）
- 每个结局通过 associations 端点分批查询 28 个 SNP

## 3. 结果（IVW）

| 结局 | nsnp | beta | SE | P | OR | 95% CI |
|---|---|---|---|---|---|---|
| CHD (CARDIoGRAM 2011) | 25 | 0.577 | 0.304 | 0.058 | 1.780 | 0.981–3.229 |
| CAD (2022) | 28 | 0.291 | 0.255 | 0.254 | 1.338 | 0.812–2.204 |
| Multiple sclerosis | 27 | -0.125 | 0.197 | 0.527 | 0.883 | 0.600–1.299 |
| Lymphoid leukaemia | 28 | -0.001 | 0.001 | 0.079 | 0.999 | 0.997–1.000 |
| HOMA-IR (T2D 相关) | 25 | -0.046 | 0.045 | 0.312 | 0.955 | 0.874–1.044 |
| Oestradiol | 28 | -0.037 | 0.046 | 0.417 | 0.964 | 0.881–1.054 |

（加权中位数：CAD OR=1.53, P=3.4e-09；CHD 2011 OR=1.94, P=6.8e-04 显著）

## 4. 解读

- LDL-C 对 CHD/CAD 呈现**一致的升高风险方向**（OR>1），
  且加权中位数法在 CHD（OR=1.94, P=6.8e-4）与 CAD（OR=1.53, P=3.4e-9）均显著
- 对多发性硬化、淋巴白血病、HOMA-IR、雌二醇未见显著效应（OR≈1）——
  体现**效应特异性**（LDL-C 主要影响心血管结局）
- 森林图见 `04.figures/mr_scan_forest.png`

## 5. 结论与教学要点

1. MR 扫描 = 单暴露 × 多结局的批量 MR，是 PheWAS 的因果推断版本
2. 显著性与方向一致性联合判读：效应应集中在生物学相关结局
3. 多结局多重比较：建议控制 FDR（实际批量研究必备）
4. 可扩展：全部结局数据集（50,000+）批量扫描 → 完整 PheWAS-MR

## 6. 可复现

```bash
# 结局关联数据（22_opengwas_api.py associations, 已入库 02.analysis/mr_scan/）
Rscript scripts/34_mr_scan.R
```
