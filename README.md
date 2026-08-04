# MR_101 — 孟德尔随机化(Mendelian Randomization)学习仓库

> 以 Git 开发规范管理的一个 MR 分析学习项目：从数据准备、工具搭建、分析流程到
> 报告产出，覆盖**两样本 MR 主流程 + 六类扩展分析**（单样本 MR、MVMR、中介 MR、
> 径向 MR、MRMix、共定位）的完整实践。

## 项目目标

- 系统学习孟德尔随机化的核心原理与三大假设（相关性、独立性、排他性）
- 掌握 TwoSampleMR / MendelianRandomization / MR-PRESSO / RadialMR / MRMix / coloc 等主流工具
- 以「端粒长度 → 冠心病 (CHD)」真实数据跑通两样本 MR 完整流程
- 以模拟+真实数据覆盖六类进阶 MR 分析类型
- 每个功能模块（工具、数据、流程、扩展）均产出可复现的测试与报告

## 目录结构

```
MR_101/
├── README.md            # 项目说明（本文件）
├── docs/                # 学习与报告文档
│   ├── 00_learning_plan.md    # 学习计划
│   ├── 01_tool_test.md        # 工具安装与测试报告
│   ├── 02_data_prep.md        # 数据准备报告
│   ├── 03_analysis.md         # 两样本 MR 分析流程与结果报告
│   ├── 04_final_report.md     # 最终学习总结报告
│   ├── 05_one_sample_mr.md    # 扩展1: 单样本 MR 报告
│   ├── 06_mvmr.md             # 扩展2: 多变量 MR 报告
│   ├── 07_mediation_mr.md     # 扩展3: 两步中介 MR 报告
│   ├── 08_radial_mr.md        # 扩展4: 径向 MR 报告
│   ├── 09_mrmix.md            # 扩展5: MRMix 报告
│   ├── 10_knowledge_base.md   # MR 知识库（全图谱 + STROBE-MR）
│   ├── 11_coloc.md            # 扩展6: 共定位分析报告
│   ├── 12_opengwas.md         # 深化: OpenGWAS 真实数据 LDL-C→CHD 报告
│   ├── 13_opengwas_guide.md   # OpenGWAS 数据库使用说明（数据内容/API/工具）
│   ├── 14_opengwas_online.md  # OpenGWAS 在线分析报告（LDL-C→CHD + 多暴露对比）
│   ├── 15_nonlinear_mr.md     # 进阶1: 非线性 MR 报告
│   ├── 16_drug_target_mr.md   # 进阶2: 药物靶点 MR（cis-MR）报告
│   ├── 17_bayesian_mr.md      # 进阶3: 贝叶斯 MR（cML）报告
│   ├── 18_mvmr_mediation.md   # 进阶4: 多变量中介 MR 报告
│   ├── 19_mr_scan.md          # 进阶5: MR 扫描（多结局批量）报告
│   ├── 20_visual_report.md    # 进阶6: MR 可视化报告（森林图/热图）
│   ├── 21_eqtl_mapping.md     # 新模块1: eQTL 定位（转录组 MR + coloc）
│   ├── 22_multiomics_mr.md    # 新模块2: 多组学 MR（eQTL/pQTL/代谢组）
│   ├── 23_drug_repurposing.md # 新模块3: 药物重定位扫描报告
│   ├── 24_bidirectional_mr.md # 补充1: 双向 MR 报告
│   ├── 25_mr_power.md         # 补充2: MR 功效计算报告
│   └── 26_pathway_mr.md       # 补充3: 通路 MR（顺序中介）报告
├── 00.data/             # 暴露/结局 GWAS 汇总数据（不入库）
├── 01.tools/            # 工具依赖与版本说明
├── 02.analysis/         # 分析中间产物与结果（CSV/绘图数据）
├── 03.reports/          # 渲染后的报告（HTML/Markdown）
├── 04.figures/          # 图表输出
└── scripts/             # 可复现分析脚本（R）
    ├── 00_env.R         # 环境检查
    ├── 01_data_prep.R   # 数据准备：工具变量提取、LD clumping、harmonise
    ├── 02_mr_analysis.R # 两样本主分析：IVW / Egger / 加权中位数等五种方法
    ├── 03_sensitivity.R # 敏感性分析：异质性、多效性、leave-one-out、PRESSO
    ├── 10_one_sample_mr.R   # 扩展1: 单样本 MR（2SLS）
    ├── 11_mvmr.R            # 扩展2: 多变量 MR（MV-IVW/Egger/Median）
    ├── 12_mediation_mr.R    # 扩展3: 两步中介 MR
    ├── 13_radial_mr.R       # 扩展4: 径向 MR 离群检验
    ├── 14_mrmix.R           # 扩展5: MRMix 稳健混合模型
    ├── 15_coloc.R           # 扩展6: 共定位分析 coloc.abf
    ├── 30_nonlinear_mr.R    # 进阶1: 非线性 MR（分位数分层）
    ├── 31_drug_target_mr.R  # 进阶2: 药物靶点 MR（cis-MR, PCSK9）
    ├── 32_bayesian_mr.R     # 进阶3: 贝叶斯 MR（cML 约束最大似然）
    ├── 33_mvmr_mediation.R  # 进阶4: 多变量中介 MR（MVMR 校正中介）
    ├── 34_mr_scan.R         # 进阶5: MR 扫描（LDL-C × 6 结局批量）
    ├── 35_visual_report.R   # 进阶6: MR 可视化报告（森林图/热图）
    ├── 36_eqtl_mapping.R    # 新模块1: eQTL 定位（HMGCR 表达→CHD + coloc）
    ├── 37_multiomics_mr.R   # 新模块2: 多组学 MR（eQTL/pQTL/代谢组→CHD）
    ├── 38_drug_repurposing.R # 新模块3: 药物重定位扫描（PCSK9/HMGCR/CETP/NPC1L1→CHD）
    ├── 39_bidirectional_mr.R # 补充1: 双向 MR（LDL-C↔CHD 反向因果）
    ├── 40_mr_power.R         # 补充2: MR 功效计算（功率曲线与样本量）
    └── 41_pathway_mr.R       # 补充3: 通路 MR（X→M1→M2→Y 顺序中介）
```

## 快速开始

```bash
# 1. 环境检查
Rscript scripts/00_env.R

# 2. 数据准备（内置真实示例数据 telomere_length -> chd）
Rscript scripts/01_data_prep.R

# 3. 主分析 + 敏感性分析
Rscript scripts/02_mr_analysis.R
Rscript scripts/03_sensitivity.R

# 4. 扩展分析（六类进阶 MR）
Rscript scripts/10_one_sample_mr.R   # 单样本 MR
Rscript scripts/11_mvmr.R            # 多变量 MR
Rscript scripts/12_mediation_mr.R    # 两步中介 MR
Rscript scripts/13_radial_mr.R       # 径向 MR
Rscript scripts/14_mrmix.R           # MRMix
Rscript scripts/15_coloc.R           # 共定位

# 5. OpenGWAS 真实数据深化（LDL-C -> CHD）
Rscript scripts/20_opengwas_harmonise.R   # 本地下载版: 工具变量提取 + harmonise
Rscript scripts/21_opengwas_analysis.R    # 本地下载版: 五方法 + 敏感性分析

# 6. OpenGWAS 在线 API 分析（需 token, 见 docs/13_opengwas_guide.md）
python3 scripts/22_opengwas_api.py gwasinfo --ids ieu-b-110,ieu-a-7   # 元数据
python3 scripts/22_opengwas_api.py tophits --id ieu-b-110 --p1 5e-8    # 工具变量
python3 scripts/22_opengwas_api.py associations --id ieu-a-7 --snps ... # 结局关联
python3 scripts/22_opengwas_api.py phewas --snps rs11591147            # PheWAS 扫描
Rscript scripts/23_opengwas_online.R   # 在线版 LDL-C -> CHD 全面分析
Rscript scripts/24_opengwas_multi.R    # 在线版 LDL/HDL/TG -> CHD 多暴露对比

# 7. 进阶分析（知识库进阶方向）
Rscript scripts/30_nonlinear_mr.R      # 非线性 MR（分位数分层）
Rscript scripts/31_drug_target_mr.R    # 药物靶点 MR（cis-MR, PCSK9）
Rscript scripts/32_bayesian_mr.R       # 贝叶斯 MR（cML）
Rscript scripts/33_mvmr_mediation.R    # 多变量中介 MR
Rscript scripts/34_mr_scan.R           # MR 扫描（LDL-C × 6 结局批量）
Rscript scripts/35_visual_report.R     # MR 可视化报告（森林图/热图）

# 8. 新模块（多组学与药物研发）
Rscript scripts/36_eqtl_mapping.R      # eQTL 定位（HMGCR 表达→CHD + coloc）
Rscript scripts/37_multiomics_mr.R     # 多组学 MR（eQTL/pQTL/代谢组→CHD）
Rscript scripts/38_drug_repurposing.R  # 药物重定位扫描（4 靶点→CHD）

# 9. 补充分析（双向 MR/功效计算/通路 MR）
Rscript scripts/39_bidirectional_mr.R  # 双向 MR（LDL-C↔CHD 反向因果）
Rscript scripts/40_mr_power.R          # MR 功效计算（功率曲线与样本量）
Rscript scripts/41_pathway_mr.R        # 通路 MR（X→M1→M2→Y 顺序中介）
```

## 结果速览

- **两样本 MR**：端粒长度 → 冠心病，17 个工具变量；
  加权中位数 OR=0.70 (P=0.023)，敏感性分析全部通过（`02.analysis/`、`docs/03_analysis.md`）
- **扩展 1 单样本 MR**：2SLS 效应 0.583（真实 0.5），OLS 有偏 0.856（`docs/05`）
- **扩展 2 MVMR**：MV-IVW X1=0.371 / X2=-0.272（真实 0.4/-0.3）（`docs/06`）
- **扩展 3 中介 MR**：间接效应 0.31，中介比例 49.2%（`docs/07`）
- **扩展 4 径向 MR**：检出离群 rs2736428，剔除后 P=6e-05（`docs/08`）
- **扩展 5 MRMix**：θ=0.14 (P=1e-17)，π0=0.46（`docs/09`）
- **扩展 6 共定位**：三区域 PP.H4>0.99（`docs/11`）
- **OpenGWAS 深化**：真实数据 LDL-C → CHD，64 个工具变量，
  IVW OR=1.57 (P=6.6e-20)，五方法全显著、Egger 截距 P=0.39 无多效性（`docs/12`）
- **OpenGWAS 在线**：API 在线版 LDL-C→CHD（加权中位数 OR=1.77, P=1.4e-06）；
  多暴露对比 LDL OR=1.77 / HDL OR=0.66（保护）/ TG OR=1.21（`02.analysis/opengwas/online*`）
- **进阶1 非线性 MR**：分位数分层揭示层间效应差异（`docs/15`）
- **进阶2 药物靶点 MR**：PCSK9→LDL-C P=8e-20（cis-MR, `docs/16`）
- **进阶3 贝叶斯 MR**：cML beta=0.47, P=0.0017 比 IVW 更稳健（`docs/17`）
- **进阶4 多变量中介 MR**：中介比例 53.6%（真实 50%, `docs/18`）
- **进阶5 MR 扫描**：LDL-C × 6 结局，CAD 加权中位数 OR=1.53 (P=3.4e-09)（`docs/19`）
- **进阶6 可视化报告**：森林图 + 多暴露/多结局热图（`docs/20`）
- **新模块1 eQTL 定位**：HMGCR 表达→CHD 转录组 MR + coloc（`docs/21`）
- **新模块2 多组学 MR**：代谢组 LDL-C OR=1.74 (P=3e-56) 显著（`docs/22`）
- **新模块3 药物重定位**：PCSK9/HMGCR/CETP/NPC1L1 四靶点→CHD 扫描（`docs/23`）
- **补充1 双向 MR**：LDL-C→CHD（OR=1.50, P=0.067）vs CHD→LDL-C（OR=1.12, P=2.3e-150）（`docs/24`）
- **补充2 功效计算**：LDL-C→CHD 实际功效仅 10.8%（`docs/25`）
- **补充3 通路 MR**：X→M1→M2→Y 双中介串联路径估计 0.110（真实 0.12）（`docs/26`）
- 说明：OpenGWAS API 需 JWT token（见 docs/13_opengwas_guide.md），本仓库
  Python 客户端（22_opengwas_api.py）负责数据获取、R 脚本（23/24）负责分析

## Git 开发规范

- 分支：`main` 为主干，功能开发用特性分支后合并
- 提交：遵循 Conventional Commits（`feat:` / `docs:` / `test:` / `fix:`）
- 报告与文档随代码同步提交，保证全流程可复现
