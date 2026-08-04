# MR_101 — 孟德尔随机化(Mendelian Randomization)学习仓库

> 以 Git 开发规范管理的一个 MR 分析入门项目：从数据准备、工具搭建、分析流程到
> 报告产出，覆盖两样本孟德尔随机化（Two-sample MR）的完整实践。

## 项目目标

- 系统学习孟德尔随机化的核心原理与三大假设（相关性、独立性、排他性）
- 掌握 TwoSampleMR / MendelianRandomization / MR-PRESSO 等主流工具的用法
- 以「LDL 胆固醇 → 冠心病 (CHD)」经典实例跑通完整分析流程
- 每个功能模块（工具、数据、流程）均产出可复现的测试与报告

## 目录结构

```
MR_101/
├── README.md            # 项目说明（本文件）
├── docs/                # 学习与报告文档
│   ├── 00_learning_plan.md    # 学习计划
│   ├── 01_tool_test.md        # 工具安装与测试报告
│   ├── 02_data_prep.md        # 数据准备报告
│   ├── 03_analysis.md         # 分析流程与结果报告
│   └── 04_final_report.md     # 最终学习总结报告
├── 00.data/             # 暴露/结局 GWAS 汇总数据（不入库）
├── 01.tools/            # 工具依赖与版本说明
├── 02.analysis/         # 分析中间产物与结果（CSV/绘图数据）
├── 03.reports/          # 渲染后的报告（HTML/Markdown）
├── 04.figures/          # 图表输出
└── scripts/             # 可复现分析脚本（R）
    ├── 00_env.R         # 环境检查
    ├── 01_data_prep.R   # 数据准备：工具变量提取、LD clumping、harmonise
    ├── 02_mr_analysis.R # 主分析：IVW / Egger / 加权中位数等五种方法
    └── 03_sensitivity.R # 敏感性分析：异质性、多效性、leave-one-out、PRESSO
```

## 快速开始

```bash
# 1. 环境检查
Rscript scripts/00_env.R

# 2. 数据准备（TwoSampleMR 内置示例数据 ldlc -> chd）
Rscript scripts/01_data_prep.R

# 3. 主分析 + 敏感性分析
Rscript scripts/02_mr_analysis.R
Rscript scripts/03_sensitivity.R
```

## 结果速览

- 暴露：LDL 胆固醇 (LDL-C)，工具变量数见 `02.analysis/instruments.txt`
- 结局：冠心病 (CHD)
- 主要方法结果：`02.analysis/mr_results.csv`
- 敏感性分析：`02.analysis/sensitivity/`
- 图形：`04.figures/`

## Git 开发规范

- 分支：`main` 为主干，功能开发用特性分支后合并
- 提交：遵循 Conventional Commits（`feat:` / `docs:` / `test:` / `fix:`）
- 报告与文档随代码同步提交，保证全流程可复现
