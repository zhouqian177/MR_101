#!/usr/bin/env Rscript
# 00_smoke_test.R — 工具链冒烟测试
# 用 TwoSampleMR 包自带的真实 GWAS 汇总数据（GIANT BMI 暴露 + CARDIoGRAM CHD 结局）
# 跑通 格式化 -> harmonise -> 五种方法 MR -> 敏感性分析 -> 绘图 全流程，
# 验证工具链可用性。输出见 01.tools/smoke_test_output.txt 与 04.figures/。

options(width = 150)
extdata <- system.file("extdata", package = "TwoSampleMR")
out_file <- "01.tools/smoke_test_output.txt"
dir.create("01.tools", showWarnings = FALSE)
dir.create("04.figures", showWarnings = FALSE)
sink(out_file, type = "output", split = TRUE)

cat("========================================\n")
cat("TwoSampleMR 工具链冒烟测试\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("TwoSampleMR 版本:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("========================================\n\n")

suppressMessages(library(TwoSampleMR))

# ---------- 1. 读取并格式化暴露/结局数据 ----------
cat("[1] 读取内置真实数据\n")
# 暴露: telomere_length.txt（端粒长度 GWAS, Codd et al. 2013）
# 结局: cardiogram.txt （冠心病 CHD, CARDIoGRAM 联盟）
# 两文件为 TwoSampleMR 官方示例数据，31 个 SNP 完全匹配
exposure_raw <- read.table(file.path(extdata, "telomere_length.txt"), header = TRUE, stringsAsFactors = FALSE)
outcome_raw  <- read.table(file.path(extdata, "cardiogram.txt"), header = TRUE, stringsAsFactors = FALSE)
cat("    暴露(telomere_length.txt) SNP 数:", nrow(exposure_raw),
    "| 结局(cardiogram.txt) SNP 数:", nrow(outcome_raw), "\n")

# 暴露文件无 pval 列，由 beta/se 推算（标准做法: p = 2*pnorm(-|z|)）
exposure_raw$pval <- 2 * pnorm(-abs(exposure_raw$beta / exposure_raw$se))

expo <- format_data(exposure_raw, type = "exposure", header = TRUE, phenotype_col = "Phenotype",
                    snp_col = "SNP", beta_col = "beta", se_col = "se",
                    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                    eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize")
outc <- format_data(outcome_raw, type = "outcome", header = TRUE, phenotype_col = "Phenotype",
                    snp_col = "SNP", beta_col = "beta", se_col = "se",
                    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                    eaf_col = "eaf")
cat("    格式化后: 暴露工具变量", nrow(expo), "个, 结局匹配 SNP", nrow(outc), "个\n\n")

# ---------- 2. harmonise ----------
cat("[2] harmonise（统一效应等位基因/方向）\n")
dat <- harmonise_data(expo, outc)
cat("    harmonise 后 SNP 数:", nrow(dat), "\n")
cat("    mr_keep=TRUE 的 SNP 数:", sum(dat$mr_keep), "\n\n")

# ---------- 3. 五种方法 MR ----------
cat("[3] 五种 MR 方法\n")
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                               "mr_weighted_mode", "mr_simple_mode"))
print(res[, c("method", "nsnp", "b", "se", "pval")], row.names = FALSE)
cat("\n")

# ---------- 4. 敏感性分析 ----------
cat("[4] 敏感性分析\n")
het <- mr_heterogeneity(dat)
cat("    异质性 Cochran Q: Q =", round(het$Q[1], 2), ", P =", format(het$Q_pval[1], digits = 4), "\n")
plt <- mr_pleiotropy_test(dat)
cat("    Egger 截距 =", round(plt$egger_intercept[1], 5),
    ", P =", format(plt$pval[1], digits = 4), "(P<0.05 提示水平多效性)\n")
ss  <- mr_singlesnp(dat)
loo <- mr_leaveoneout(dat)
cat("    单 SNP 分析:", nrow(ss), "条; leave-one-out:", nrow(loo), "条\n")
write.csv(res, "02.analysis/mr_results_smoke.csv", row.names = FALSE)

# ---------- 5. 绘图 ----------
cat("[5] 绘图（散点/森林/漏斗）\n")
p1 <- mr_scatter_plot(res, dat)
p2 <- mr_forest_plot(mr_singlesnp(dat))
p3 <- mr_funnel_plot(mr_singlesnp(dat))
ggplot2::ggsave("04.figures/smoke_scatter.png", p1[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/smoke_forest.png", p2[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/smoke_funnel.png", p3[[1]], width = 7, height = 6)
cat("    已保存 04.figures/smoke_{scatter,forest,funnel}.png\n\n")

cat("========================================\n")
cat("冒烟测试完成：工具链可用 ✔\n")
sink()
cat("输出已写入:", out_file, "\n")
