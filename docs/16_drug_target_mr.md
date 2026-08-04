# 药物靶点 MR 报告（16_drug_target_mr）

> 对应模块：进阶 2 —— 药物靶点孟德尔随机化（Drug-target MR / cis-MR）。
> 脚本：`scripts/31_drug_target_mr.R`；产物：`02.analysis/drug_target/`

## 1. 方法学

药物靶点 MR 用**靶点蛋白/基因表达**的 pQTL/eQTL 变异作为工具变量，评估
"调节该靶点"对结局的因果效应，从而**模拟药物疗效**（药物重定位、靶点验证）。

- **cis-MR**：只用靶点基因 ±1Mb 内的顺式变异（cis-pQTL/cis-eQTL），
  多效性风险更低、因果方向更可靠
- **Wald ratio**：单 SNP 时 β_因果 = β_结局 / β_暴露
- 对二分类结局换算 OR = exp(β)，解读为"每升高 1 单位蛋白水平的风险比"

## 2. 实例：PCSK9（经典降脂药靶点）

PCSK9 是依洛尤单抗（Evolocumab）、阿利西尤单抗（Alirocumab）等
PCSK9 抑制剂的靶点。工具变量 rs631220（chr1:55527479，**PCSK9 基因内**，
cis-pQTL P=4.8e-12，数据来自 OpenGWAS ebi-a-GCST90010246）。

## 3. 结果（Wald ratio）

| 关联 | beta | SE | P | OR(每单位蛋白) | 95% CI |
|---|---|---|---|---|---|
| PCSK9 → LDL-C | 0.0625 | 0.0069 | **8.1e-20** | — | — |
| PCSK9 → CHD | 0.0134 | 0.0190 | 0.48 | 1.014 | 0.976–1.052 |

**解读**：
- **PCSK9 蛋白水平每升高 1 单位，LDL-C 升高 0.0625（P=8e-20）**——
  支持 PCSK9 靶点对 LDL-C 的强因果效应
- PCSK9 → CHD 方向为正（升高 PCSK9 增加 CHD 风险）但不显著（P=0.48），
  符合预期：PCSK9 对 CHD 的效应主要通过 LDL-C 介导，单 SNP 功效有限
- 临床意义：抑制 PCSK9 → 降低 LDL-C → 降低 CHD 风险（与临床试验一致）

## 4. 结论与教学要点

1. 药物靶点 MR 是**药物研发早期证据**的重要手段（靶点验证/重定位）
2. cis 变异多效性风险低，但单变异工具功效有限 → 可扩展多 cis-pQTL 联合
3. 常见应用：pQTL（蛋白）、eQTL（基因表达）、sQTL 等作为靶点工具
4. 数据：OpenGWAS 有大量 pQTL 数据集（如 prot-a-*、ebi-a-GCST90* 系列），
   可用 `python3 scripts/22_opengwas_api.py search --trait PCSK9` 检索

## 5. 可复现

```bash
# 在线获取 PCSK9 pQTL 与结局关联（见 docs/13_opengwas_guide.md）
Rscript scripts/31_drug_target_mr.R
```
