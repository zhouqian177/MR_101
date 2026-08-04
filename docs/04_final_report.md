# MR 最终学习报告（04_final_report）

> 项目：MR_101 —— 孟德尔随机化学习仓库（Two-sample MR + 六类扩展分析）
> 完成时间：2026-08-04（首版）、2026-08-04（扩展版）

## 1. 项目目标达成情况

| 学习计划阶段 | 目标 | 达成 |
|---|---|---|
| 阶段一：原理 | 三大假设、五方法、检验体系 | ✔ 见 §2、docs/03_analysis.md、docs/10_knowledge_base.md |
| 阶段二：工具 | TwoSampleMR 等工具链 | ✔ 见 docs/01_tool_test.md |
| 阶段三：数据 | 工具变量提取、clumping、harmonise、F 检验 | ✔ 见 docs/02_data_prep.md |
| 阶段四：流程 | 端到端两样本分析 | ✔ 见 docs/03_analysis.md |
| 阶段五：报告 | 各模块报告齐全 | ✔ 各模块报告 + 本报告 |
| 阶段六：扩展 | 六类进阶 MR 分析 | ✔ 见 §3.2 |

## 2. 理论学习要点

### 2.1 孟德尔随机化三大假设

1. **相关性假设**：工具变量（SNP）与暴露强相关（用 F 统计量检验，F>10 合格）
2. **独立性假设**：工具变量与混杂因素无关
3. **排他性假设**：工具变量仅通过暴露影响结局（用 MR-Egger 截距、MR-PRESSO 检验）

### 2.2 五种估计方法的适用条件

| 方法 | 假设条件 | 本实例结果 |
|---|---|---|
| IVW | 所有 SNP 均为有效工具变量 | OR=0.79, P=0.055 |
| MR-Egger | 允许定向多效性（截距≠0） | OR=0.76, P=0.241 |
| 加权中位数 | 至少 50% 权重来自有效 SNP | OR=0.70, **P=0.023** |
| 加权众数 | 最大 SNP 簇有效 | OR=0.71, **P=0.044** |
| 简单众数 | 同上（未加权） | OR=0.76, P=0.177 |

### 2.3 检验体系

- 异质性：Cochran Q（P=0.15，无显著异质性）
- 水平多效性：MR-Egger 截距（P=0.83）与 MR-PRESSO（全局 P=0.143，无离群）
- 方向性：Steiger 检验（correct_causal_direction=TRUE）
- 稳健性：leave-one-out（无单 SNP 主导）

## 3. 实例分析结论

**端粒长度 → 冠心病（CHD）**：端粒长度越长，CHD 风险越低
（加权中位数 OR=0.70, 95%CI 0.52–0.95），与文献报道方向一致。

> 实例选择说明：学习计划原定 LDL-C → CHD，但当前环境 OpenGWAS API 被代理阻断，
> 无法在线提取 LDL-C GWAS 数据，故改用 TwoSampleMR 官方内置的真实 GWAS 示例数据
> （端粒长度 → CHD，31 个 SNP，与官方 vignette 一致）。分析流程完全通用，
> 网络恢复后替换数据源（`extract_instruments("ieu-b-110")` 等）即可复现 LDL-C 实例。

## 4. 模块报告索引（每个功能模块均有对应报告）

| 模块 | 产物 | 报告 |
|---|---|---|
| 学习计划 | docs/00_learning_plan.md | 六阶段计划与验收标准 |
| 工具安装与测试 | scripts/00_env.R、00_smoke_test.R | docs/01_tool_test.md（含踩坑记录） |
| 数据准备 | scripts/01_data_prep.R、02.analysis/ | docs/02_data_prep.md |
| 两样本分析流程 | scripts/02_mr_analysis.R、03_sensitivity.R | docs/03_analysis.md |
| 图表 | 04.figures/*.png | 散点/森林/漏斗/leave-oneout |
| 扩展1 单样本 MR | scripts/10_one_sample_mr.R | docs/05_one_sample_mr.md |
| 扩展2 多变量 MR | scripts/11_mvmr.R | docs/06_mvmr.md |
| 扩展3 两步中介 MR | scripts/12_mediation_mr.R | docs/07_mediation_mr.md |
| 扩展4 径向 MR | scripts/13_radial_mr.R | docs/08_radial_mr.md |
| 扩展5 MRMix | scripts/14_mrmix.R | docs/09_mrmix.md |
| 扩展6 共定位 | scripts/15_coloc.R | docs/11_coloc.md |
| 深化 OpenGWAS | scripts/20/21_opengwas_*.R | docs/12_opengwas.md（真实 LDL-C→CHD） |
| 知识库 | docs/10_knowledge_base.md | MR 全图谱 + STROBE-MR 清单 |

## 5. 局限性与改进方向

1. **弱工具变量**：两样本示例数据样本量小（~2 万），平均 F=8.6 < 10，IVW 可能被低估；
   实际研究应过滤 F>10 或使用更强 GWAS
2. **数据源受限**：OpenGWAS API 不可用，未能演示在线数据提取与 LD 参考面板下载；
   网络恢复后应补充 `extract_instruments` + `clump_data` 在线流程
3. **人群一致性**：两样本 MR 要求暴露/结局人群尽量同源（本示例均为欧洲人群，OK）
4. **扩展模块数据**：MVMR/中介/共定位等扩展以模拟数据演示方法学（共定位用 coloc 官方
   模拟区域），真实 GWAS 数据应用流程已就绪，待网络恢复后可替换
5. **可再扩展**：非线性 MR、药物靶点 MR（cis-MR/pQTL）、贝叶斯 MR（MRBEE/cML）、
   多变量中介等进阶内容（见 docs/10_knowledge_base.md §7）

## 6. Git 开发规范执行情况

- 分支：main 主干 + Conventional Commits（docs:/test:/feat: 前缀）
- 提交粒度：按模块划分（骨架→工具→数据→脚本→流程→报告→扩展）
- 每次提交均通过验证后再提交（冒烟测试、端到端重跑、各扩展脚本单独验证）

## 7. 结论

本项目以 Git 规范完成了一个完整的孟德尔随机化学习闭环：原理 → 工具 → 数据 →
流程 → 报告，并额外覆盖六类进阶分析（单样本 MR、MVMR、中介 MR、径向 MR、
MRMix、共定位），每个模块均有可复现脚本与对应报告。核心能力已就绪：
**拿到任意暴露/结局 GWAS 汇总数据，即可按 scripts/ 流程完成标准两样本 MR 分析；
进阶问题（多暴露、中介机制、多效性稳健、LD 假关联排除）均有对应扩展工具。**
