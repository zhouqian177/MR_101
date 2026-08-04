# 00.data — 数据目录说明

存放暴露与结局的 GWAS 汇总数据（sumstats）。大型数据不入库（见 `.gitignore`），
通过 `scripts/01_data_prep.R` 可再生。

数据来源：

| 数据集 | 说明 | 获取方式 |
|---|---|---|
| ldlc | LDL 胆固醇 GWAS（暴露） | `TwoSampleMR::ldlc` 内置 / OpenGWAS `ieu-b-110` |
| hdlc | HDL 胆固醇 GWAS | `TwoSampleMR::hdlc` 内置 |
| tg | 甘油三酯 GWAS | `TwoSampleMR::tg` 内置 |
| chd | 冠心病 GWAS（结局） | `TwoSampleMR::chd` 内置 / OpenGWAS `ieu-a-7` |

OpenGWAS API: https://gwas.mrcieu.ac.uk/ （`ieugwasr` 包访问）
