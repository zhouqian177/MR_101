# OpenGWAS 数据库使用说明（13_opengwas_guide）

> 面向本仓库使用者的 OpenGWAS（IEU OpenGWAS）数据库完整使用指南：
> 数据库内容、可下载数据、API 访问方式、R/Python 工具与本仓库客户端。
> 更新日期：2026-08-04

## 1. 数据库简介

**IEU OpenGWAS**（https://gwas.mrcieu.ac.uk/，新门户 https://opengwas.io/）
由英国布里斯托大学 MRC IEU 维护，是孟德尔随机化研究最常用的 GWAS 汇总数据库：

- **规模**：50,000+ 个经 QC 与标准化整理的完整 GWAS 汇总数据集
- **数据来源**：公共 GWAS 文献（GLGC、CARDIoGRAMplusC4D、UKB、FinnGen、GIANT 等）
- **覆盖性状**：疾病（冠心病、2 型糖尿病、癌症…）、血液指标（血脂、血糖…）、
  人体测量（BMI、身高）、生活方式（吸烟、饮酒）、多组学（eQTL、pQTL、代谢物）等
- **人群**：主要为欧洲人群（European），部分亚洲/混合人群
- **坐标版本**：大多为 GRCh37（build 37），查询时注意

## 2. 能下载哪些数据

### 2.1 数据集元数据（gwasinfo）

每个数据集有唯一 ID（如 `ieu-b-110`），可查询：
样本量（ncase/ncontrol）、SNP 数、性状名、人群、构建版本、发表信息等。

### 2.2 数据集内汇总统计（关联数据）

- **工具变量（tophits）**：某数据集中 P 小于阈值的显著 SNP（可带 LD clumping）
- **指定 SNP 关联（associations）**：任意 SNP 在某数据集中的效应量
  （beta/se/p/eaf/效应等位基因等）
- **区域查询（variants）**：某染色体区间的全部 SNP 关联（用于 fine-mapping/共定位）
- **全量数据文件**：VCF 格式完整汇总数据（网页端生成 2 小时有效下载链接）

### 2.3 跨数据集数据

- **PheWAS**：单个变异在所有数据集中的关联扫描
- **LD 信息**：SNP 间 LD 矩阵、LD proxy（服务端或本地）
- **dbSNP 注释**：rsID ↔ 染色体位置转换、基因注释

## 3. 访问方式总览

| 方式 | 适用场景 | 认证 |
|---|---|---|
| 网页 https://api.opengwas.io/ | 人工浏览、单次查询 | 注册账号 |
| REST API | 程序化批量访问 | JWT token |
| R 包 ieugwasr / TwoSampleMR | R 工作流 | JWT token |
| Python 包 ieugwaspy | Python 工作流 | JWT token |
| 本仓库客户端 scripts/22_opengwas_api.py | 本仓库标准通道 | JWT token（.env） |

## 4. 认证：JWT Token

**2024-05-01 起所有数据端点强制要求 JWT token**（仅 /api/status 免认证）。

1. 在 https://api.opengwas.io/ 用 GitHub 或邮箱注册账号
2. 登录后在 Account 页生成 token（**有效期 14 天**，到期重新生成）
3. 配置方式（二选一）：
   ```bash
   # 方式 A: 环境变量
   export OPENGWAS_JWT="eyJhbGciOi..."
   # 方式 B: 项目 .env 文件（本仓库约定, 已被 .gitignore 排除）
   echo 'OpenGWAS_JWT="eyJhbGciOi..."' > .env
   ```

> 安全提醒：token 等同密码，**不要**提交到 git、不要粘贴到聊天中；
> 本仓库 .env 已在 .gitignore 中排除。

### 额度（Allowance）

- 注册后为 Trial/Standard 档：**100,000 次/10 分钟**（免费）
- 超额返回 429，需等待 Retry-After；频繁超限可能封 IP
- 客户端已内置 429 自动等待逻辑

## 5. 常用数据 ID（本仓库已验证）

| ID | 性状 | 人群 | 样本量 | 用途 |
|---|---|---|---|---|
| ieu-b-110 | LDL 胆固醇 | 欧洲 | 188,578 | 暴露（本仓库 LDL→CHD 分析） |
| ieu-b-109 | HDL 胆固醇 | 欧洲 | 188,578 | 暴露（多暴露对比） |
| ieu-b-111 | 甘油三酯 TG | 欧洲 | 188,578 | 暴露（多暴露对比） |
| ieu-a-7 | 冠心病 CHD | 混合 | 60,801 病例/123,504 对照 | 结局 |

> 查找更多数据集：网页端搜索，或
> `python3 scripts/22_opengwas_api.py search --trait "coronary"`。

## 6. 本仓库 Python 客户端用法

```bash
# 元数据
python3 scripts/22_opengwas_api.py gwasinfo --ids ieu-b-110,ieu-a-7

# 按关键词搜索数据集
python3 scripts/22_opengwas_api.py search --trait LDL --out 02.analysis/opengwas/dataset_search.csv

# 工具变量提取（P<5e-8, 自动保存 CSV）
python3 scripts/22_opengwas_api.py tophits --id ieu-b-110 --p1 5e-8 \
    --out 02.analysis/opengwas/online_ieu-b-110_instruments.csv

# 指定 SNP 在某数据集的关联（自动分批, N(id)×N(snp)<=64）
python3 scripts/22_opengwas_api.py associations --id ieu-a-7 \
    --snps rs1800588,rs599839 --out 02.analysis/opengwas/outcome.csv

# PheWAS: 变异跨全部数据集扫描
python3 scripts/22_opengwas_api.py phewas --snps rs11591147 --p1 5e-8 \
    --out 02.analysis/opengwas/phewas.csv
```

## 7. R 包（ieugwasr / TwoSampleMR）用法

R 官方工作流（需先配置 `OPENGWAS_JWT` 环境变量）：

```r
library(TwoSampleMR)
# 提取工具变量（自动 LD clumping, P<5e-8, r2<0.001, 10Mb）
expo <- extract_instruments(outcomes = "ieu-b-110", clump = TRUE)
# 结局关联
outc <- extract_outcome_data(snps = expo$SNP, outcomes = "ieu-a-7")
# harmonise + MR
dat <- harmonise_data(expo, outc)
mr(dat)
```

> 本环境限制：系统 R 的 libcurl 7.61.1 过旧，无法通过代理完成 TLS 握手
> （`SSL_ERROR_SYSCALL`），故 R 包直连不可用；本仓库以 Python 客户端
> 完成数据获取、以 R 完成分析（scripts/23、24），等价且已验证。

## 8. 数据字段说明（associations/tophits 输出）

| 字段 | 含义 |
|---|---|
| rsid | 变异 rsID |
| chr / position | 染色体与位置（build37） |
| ea / nea | 效应等位基因 / 非效应等位基因 |
| eaf | 效应等位基因频率 |
| beta | 效应量（连续性状）或 log(OR)（二分类） |
| se | 标准误 |
| p | P 值 |
| n / ncase / ncontrol | 样本量 / 病例 / 对照 |

## 9. 常见问题

| 问题 | 处理 |
|---|---|
| 401 Unauthorized | token 缺失/过期 → 重新生成并配置 .env |
| 400 "N(id)*N(variant)<=64" | 分批查询（客户端已自动处理） |
| 429 Too Many Requests | 等待 Retry-After（客户端自动处理） |
| gwasinfo 返回全量 | 该端点忽略 id 参数，客户端已本地过滤 |
| tophits GET 405 | 该端点需 POST（客户端已实现） |
| 数据是 build37 坐标 | 与参考面板/结局数据需同版本 |

## 10. 与公开下载数据的对应关系

OpenGWAS 中多数数据集源自公开 GWAS，可在原机构直接下载（内容一致）：

| OpenGWAS ID | 原始来源 | 本仓库本地文件 |
|---|---|---|
| ieu-b-110 | GLGC 2013 LDL | 00.data/gwas/LDL_GLGC2013.txt.gz |
| ieu-a-7 | CARDIoGRAMplusC4D 2015 | 00.data/gwas/CHD_CAD2015.txt |

本仓库同时支持两种路径（在线 API 与本地文件），分析脚本（scripts/23、24）可直接复用。
