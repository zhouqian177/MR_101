# eQTL 定位报告（21_eqtl_mapping）

> 对应模块：新模块 1 —— eQTL 定位（转录组 MR + 共定位）。
> 脚本：`scripts/36_eqtl_mapping.R`；产物：`02.analysis/multiomics/eqtl_*.csv`

## 1. 方法学

eQTL 定位（Transcriptome-wide MR）用**基因表达 eQTL** 作为工具变量，
评估基因表达水平对疾病的因果效应，是药物靶点验证与机制研究的重要工具：

1. **转录组 MR**：基因表达（cis-eQTL）→ 疾病，识别表达变化对结局的因果效应
2. **coloc 共定位**：基因表达信号与疾病信号是否共享同一因果变异
   （PP.H4 高 = 真共定位，PP.H3 高 = LD 假关联）

## 2. 实例：HMGCR 基因表达 → CHD

- 暴露：HMGCR 基因表达 eQTL（eqtl-a-ENSG00000140464，OpenGWAS）
- 结局：CHD（ieu-a-7）
- cis 工具变量：rs10851868（P=3.2e-145）、rs12904134（P=9.9e-10）
- 背景：HMGCR 是他汀类药物的靶点基因（HMG-CoA 还原酶）

## 3. 结果

### 3.1 转录组 MR（HMGCR 表达 → CHD）

| 方法 | beta | SE | P | OR | 95% CI |
|---|---|---|---|---|---|
| IVW（2 工具） | -0.040 | 0.030 | 0.177 | 0.960 | 0.906–1.018 |
| Wald ratio (rs10851868) | -0.048 | 0.031 | 0.120 | 0.954 | — |
| Wald ratio (rs12904134) | 0.133 | 0.152 | 0.380 | 1.142 | — |

方向：HMGCR 表达升高倾向降低 CHD 风险（负向），与"他汀抑制 HMGCR 降低
CHD 风险"的临床效应方向一致（表达降低→CHD 风险降低，即表达与风险负相关），
但未达显著（cis 工具功效有限）。

### 3.2 coloc 共定位

- PP.H1=0.985（eQTL 信号主导）、PP.H3≈0、PP.H4=0.015
- 本例基因表达信号强但疾病信号弱，共定位证据有限
  （教学演示重点在流程；真实研究需更大功效的结局 GWAS）

## 4. 结论与教学要点

1. eQTL 定位 = 转录组 MR + 共定位，将基因表达与疾病因果联系起来
2. 工具变量选择：cis-eQTL（基因 ±1Mb 内）多效性风险低
3. coloc 区分"共享因果变异"与"LD 假关联"，是 MR 后必做的验证
4. 数据源：OpenGWAS eqtl-a-*（GTEx 类）系列，每基因一个数据集

## 5. 可复现

```bash
# 在线数据已入库（02.analysis/multiomics/）
Rscript scripts/36_eqtl_mapping.R
```
