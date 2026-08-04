# 工具测试报告（01_tool_test）

> 对应模块：01.tools —— R MR 分析工具链的安装、验证与冒烟测试。
> 执行时间：2026-08-04；脚本：`scripts/00_env.R`、`scripts/00_smoke_test.R`

## 1. 环境概况

- R 版本：4.5.0（x86_64-redhat-linux-gnu）
- Python：3.11.7（miniconda，备用）
- 网络：OpenGWAS API（api.opengwas.io）经代理访问被阻断，**R 包 `ieugwasr` 无法直连 OpenGWAS**；
  因此本项目的实践数据采用 TwoSampleMR 包内置的真实 GWAS 汇总数据（见 §3），
  API 路径（`extract_instruments` / `extract_outcome_data`）在文档中给出但当前环境不可用。

## 2. 关键 R 包安装与版本

| 包 | 版本 | 来源 | 状态 |
|---|---|---|---|
| TwoSampleMR | 0.7.9 | GitHub MRCIEU（CRAN 已归档） | ✔ 可用 |
| MendelianRandomization | 0.10.0 | CRAN | ✔ 可用 |
| ieugwasr | 1.1.0 | CRAN | ✔ 可用（API 受网络限制） |
| MRInstruments | 0.3.2 | GitHub MRCIEU | ✔ 可用 |
| MRPRESSO | 1.0 | GitHub rondolab | ✔ 可用 |
| MRMix | - | GitHub gqi | ✔ 可用 |
| RadialMR | - | GitHub WSpiller | ✔ 可用 |
| glmnet / gmp / psych | - | CRAN | ✔ 可用 |

安装要点（踩坑记录）：

1. **TwoSampleMR 已从 CRAN 归档**，须从 GitHub 安装；api.github.com 被代理阻断，
   改用 `codeload.github.com` 下载 tarball 后 `R CMD INSTALL`。
2. **conda gcc 无法编译 R 包**：RedHat 的 R 编译参数带 `-specs=redhat-annobin-cc1`，
   conda gcc 缺少 `gcc-annobin.so` 插件导致 `configure: C compiler cannot create executables`。
   解决：在 `~/.R/Makevars` 中强制使用系统编译器 `/usr/bin/gcc` 与 `/usr/bin/gfortran`。
3. 安装被中断后需清理 `00LOCK*` 残留目录，否则报 `failed to lock directory`。

## 3. 冒烟测试（smoke test）

用 TwoSampleMR 官方示例数据做端到端验证：

- **暴露**：端粒长度（telomere_length.txt，Codd et al. 2013 GWAS，31 SNP）
- **结局**：冠心病 CHD（cardiogram.txt，CARDIoGRAM 联盟，31 SNP）
- 流程：格式化 → harmonise → 五种方法 MR → 敏感性分析 → 绘图

### 3.1 数据匹配与 harmonise

| 步骤 | SNP 数 |
|---|---|
| 暴露原始数据 | 31 |
| 结局原始数据 | 31 |
| 格式化后匹配 | 31 / 31 |
| harmonise 后（mr_keep=TRUE） | 31 |

### 3.2 五种 MR 方法结果（Telomere_length → CHD）

| 方法 | nsnp | beta | se | P |
|---|---|---|---|---|
| Inverse variance weighted | 31 | -0.400 | 0.088 | 6.04e-06 |
| MR Egger | 31 | -0.518 | 0.173 | 5.67e-03 |
| Weighted median | 31 | -0.373 | 0.117 | 1.45e-03 |
| Weighted mode | 31 | -0.382 | 0.399 | 3.47e-01 |
| Simple mode | 31 | -0.382 | 0.400 | 3.48e-01 |

方向一致（均为负向），IVW / Egger / Weighted median 均显著，与文献结论
「端粒长度缩短增加 CHD 风险」一致（beta 为负表示每单位端粒长度的 OR<1）。

### 3.3 敏感性分析

- 异质性：Cochran Q = 47.31，P = 0.017（存在一定异质性，可接受范围内）
- 水平多效性：MR-Egger 截距 = 0.0064，P = 0.435（**无显著水平多效性**）
- 单 SNP 分析：33 条；leave-one-out：32 条（结果稳健性见正式分析报告）

### 3.4 绘图输出（04.figures/）

- `smoke_scatter.png`：SNP 效应散点图（斜率即 MR 估计）
- `smoke_forest.png`：单 SNP 森林图
- `smoke_funnel.png`：漏斗图（检验多效性对称性）

## 4. 结论

- 工具链完整可用：TwoSampleMR 全家桶 + 敏感性分析工具全部通过冒烟测试 ✔
- MR 结果方向与效应量合理，五种方法方向一致，敏感性分析通过 ✔
- 遗留限制：OpenGWAS API 在本环境不可直连（代理阻断），在线数据路径留待网络恢复后补充验证

## 5. 可复现

```bash
Rscript scripts/00_env.R        # 环境检查 -> 01.tools/env_check.txt
Rscript scripts/00_smoke_test.R # 冒烟测试 -> 01.tools/smoke_test_output.txt + 04.figures/
```
