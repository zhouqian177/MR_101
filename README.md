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
│   └── 12_opengwas.md         # 深化: OpenGWAS 真实数据 LDL-C→CHD 报告
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
    └── 15_coloc.R           # 扩展6: 共定位分析 coloc.abf
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
Rscript scripts/20_opengwas_harmonise.R   # 工具变量提取 + harmonise
Rscript scripts/21_opengwas_analysis.R    # 五方法 + 敏感性分析
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
- 说明：OpenGWAS API 需 JWT token（本环境未配置），改用与其同源的公开下载数据
  （GLGC 2013 LDL = ieu-b-110、CARDIoGRAMplusC4D 2015 = ieu-a-7），数据完全一致

## Git 开发规范

- 分支：`main` 为主干，功能开发用特性分支后合并
- 提交：遵循 Conventional Commits（`feat:` / `docs:` / `test:` / `fix:`）
- 报告与文档随代码同步提交，保证全流程可复现
