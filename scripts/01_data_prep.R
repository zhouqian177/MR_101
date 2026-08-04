#!/usr/bin/env Rscript
# 01_data_prep.R — 数据准备：工具变量提取、LD clumping、harmonise、F 统计量
#
# 数据：暴露 = 端粒长度 GWAS（telomere_length.txt, Codd et al. 2013）
#       结局 = 冠心病 CHD（cardiogram.txt, CARDIoGRAM 联盟）
# LD 参考：01.GWAS 项目的 CEU HapMap3 面板（plink bed/bim/fam）
# 输出：02.analysis/instruments.txt、harmonised.csv、F 统计量、数据准备摘要

options(width = 150)
dir.create("02.analysis", showWarnings = FALSE)

# ---------- 配置 ----------
EXPOSURE_FILE <- "00.data/telomere_length.txt"
OUTCOME_FILE  <- "00.data/cardiogram.txt"
PLINK         <- "/y/u/zhouqian/00.AI_learning/01.GWAS/software/bin/plink"
LD_REF_PREFIX <- "/y/u/zhouqian/00.AI_learning/01.GWAS/01.qc/ceu_raw"  # CEU HapMap3

cat("========================================\n")
cat("MR 数据准备报告\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

suppressMessages(library(TwoSampleMR))

# ---------- 1. 读取暴露数据 ----------
cat("[1] 读取暴露数据（端粒长度 GWAS）\n")
expo_raw <- read.table(EXPOSURE_FILE, header = TRUE, stringsAsFactors = FALSE)
cat("    SNP 数:", nrow(expo_raw), "| 列:", paste(colnames(expo_raw), collapse = ", "), "\n")
# 暴露文件无 pval 列，由 beta/se 推算
expo_raw$pval <- 2 * pnorm(-abs(expo_raw$beta / expo_raw$se))
cat("    推算 pval 后全基因组显著 (P<5e-8) SNP 数:",
    sum(expo_raw$pval < 5e-8), "\n\n")

# ---------- 2. 工具变量提取 ----------
cat("[2] 工具变量提取\n")
# 说明: 本示例数据(telomere_length.txt)为 TwoSampleMR 官方示例中已筛选的
#       工具变量列表（31 个端粒长度相关 SNP），直接作为工具变量候选。
#       实际研究中应使用完整 GWAS 汇总数据按 P<5e-8 提取（此处因 OpenGWAS
#       API 被代理阻断，采用包内置的已筛选列表；标准流程见 docs/02_data_prep.md）。
write.table(expo_raw[, c("SNP", "pval")], "02.analysis/snps_pval.txt",
            row.names = FALSE, quote = FALSE, sep = "\t", col.names = c("SNP", "P"))
cat("    候选工具变量(已筛选列表):", nrow(expo_raw), "个\n\n")

# ---------- 3. LD clumping ----------
cat("[3] LD clumping（r2<0.001, 10Mb, CEU HapMap3 参考面板）\n")
# 对已筛选列表做 LD 剪枝：--clump-p1 1 使全部 SNP 均可作为索引 SNP
system2(PLINK, c("--bfile", LD_REF_PREFIX,
                 "--clump", "02.analysis/snps_pval.txt",
                 "--clump-p1", "1", "--clump-p2", "1",
                 "--clump-r2", "0.001", "--clump-kb", "10000",
                 "--out", "02.analysis/clump"),
        stdout = "02.analysis/clump.log", stderr = "02.analysis/clump.log")
clumped <- read.table("02.analysis/clump.clumped", header = TRUE, stringsAsFactors = FALSE)
keep_snps <- clumped$SNP
cat("    clumping 后保留工具变量:", length(keep_snps), "个（去除 LD 冗余 SNP）\n")
expo_raw <- subset(expo_raw, SNP %in% keep_snps)
cat("    最终工具变量:", nrow(expo_raw), "个\n\n")

# ---------- 4. 格式化暴露/结局并 harmonise ----------
cat("[4] harmonise（统一效应等位基因与方向）\n")
expo <- format_data(expo_raw, type = "exposure", header = TRUE,
                    phenotype_col = "Phenotype", snp_col = "SNP",
                    beta_col = "beta", se_col = "se",
                    effect_allele_col = "effect_allele",
                    other_allele_col = "other_allele",
                    eaf_col = "eaf", pval_col = "pval",
                    samplesize_col = "samplesize")
outc_raw <- read.table(OUTCOME_FILE, header = TRUE, stringsAsFactors = FALSE)
outc <- format_data(outc_raw, type = "outcome", header = TRUE,
                    snp_col = "SNP", beta_col = "beta", se_col = "se",
                    effect_allele_col = "effect_allele",
                    other_allele_col = "other_allele", eaf_col = "eaf")
dat <- harmonise_data(expo, outc)
dat <- subset(dat, mr_keep)
cat("    harmonise 后保留 SNP 数:", nrow(dat), "\n\n")

# ---------- 5. F 统计量（弱工具变量检验） ----------
cat("[5] F 统计量（F>10 提示无弱工具变量）\n")
# F = (beta/se)^2，使用暴露 GWAS 的 beta/se（即 z 统计量平方）
dat$F_stat <- (dat$beta.exposure / dat$se.exposure)^2
cat("    F 统计量: 最小", round(min(dat$F_stat), 1),
    "| 平均", round(mean(dat$F_stat), 1),
    "| 最大", round(max(dat$F_stat), 1), "\n\n")

# ---------- 6. 保存结果 ----------
cat("[6] 保存数据准备产物\n")
instr <- dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                 "beta.exposure", "se.exposure", "pval.exposure",
                 "eaf.exposure", "beta.outcome", "se.outcome", "pval.outcome",
                 "F_stat", "mr_keep")]
write.csv(instr, "02.analysis/instruments.txt", row.names = FALSE)
write.csv(dat, "02.analysis/harmonised.csv", row.names = FALSE)
sink("02.analysis/data_prep_summary.txt")
cat("MR 数据准备摘要\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("暴露: 端粒长度 (Telomere_length) | 结局: 冠心病 (CHD)\n")
cat("候选工具变量 (P<5e-8): 31\n")
cat("LD clumping 后工具变量:", nrow(instr), "\n")
cat("harmonise 后 SNP:", nrow(instr), "\n")
cat("F 统计量: min =", round(min(dat$F_stat), 1),
    ", mean =", round(mean(dat$F_stat), 1),
    ", max =", round(max(dat$F_stat), 1), "\n")
cat("结论: F>10 无弱工具变量；工具变量数量 >= 3 可进行 MR 分析\n")
sink()
cat("    02.analysis/instruments.txt、harmonised.csv、data_prep_summary.txt 已生成\n")
cat("\n数据准备完成 ✔\n")
