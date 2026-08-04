#!/usr/bin/env Rscript
# 21_opengwas_analysis.R — 公开 GWAS 数据: LDL-C -> CHD 五方法 MR + 敏感性分析
#
# 输入: 02.analysis/opengwas/ldl_chd_harmonised.csv（65 个工具变量）
# 输出: 02.analysis/opengwas/ 下结果表与 04.figures/opengwas_*.png

options(width = 150)
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("LDL-C -> CHD 真实数据 MR 分析\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

d <- read.csv("02.analysis/opengwas/ldl_chd_harmonised.csv", stringsAsFactors = FALSE)
cat("[1] 工具变量数:", nrow(d), "\n")

# ---------- 格式化为 TwoSampleMR 对象 ----------
expo <- data.frame(SNP = d$rsid, beta = d$beta.exposure, se = d$se.exposure,
                   effect_allele = d$A1, other_allele = d$A2,
                   eaf = d$eaf.exposure, pval = d$pval.exposure,
                   samplesize = 188578)   # GLGC LDL 最大样本量
outc <- data.frame(SNP = d$rsid, beta = d$beta.outcome, se = d$se.outcome,
                   effect_allele = d$A1, other_allele = d$A2,
                   eaf = d$eaf.outcome, pval = d$pval.outcome,
                   samplesize = 184305)   # CARDIoGRAMplusC4D CAD
expo_f <- format_data(expo, type = "exposure", header = TRUE,
                      snp_col = "SNP", beta_col = "beta", se_col = "se",
                      effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele",
                      eaf_col = "eaf", pval_col = "pval")
outc_f <- format_data(outc, type = "outcome", header = TRUE,
                      snp_col = "SNP", beta_col = "beta", se_col = "se",
                      effect_allele_col = "effect_allele",
                      other_allele_col = "other_allele",
                      eaf_col = "eaf", pval_col = "pval")
dat <- harmonise_data(expo_f, outc_f)
dat <- subset(dat, mr_keep)
cat("[2] harmonise 后保留 SNP:", nrow(dat), "\n\n")

# ---------- 五种方法 ----------
cat("[3] 五种 MR 方法\n")
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                               "mr_weighted_mode", "mr_simple_mode"))
res$OR <- exp(res$b); res$OR_lci <- exp(res$b - 1.96 * res$se); res$OR_uci <- exp(res$b + 1.96 * res$se)
print(res[, c("method", "nsnp", "b", "se", "pval", "OR", "OR_lci", "OR_uci")], row.names = FALSE)
write.csv(res, "02.analysis/opengwas/mr_results.csv", row.names = FALSE)

# ---------- 异质性 / 多效性 ----------
cat("\n[4] 异质性与水平多效性\n")
het <- mr_heterogeneity(dat)
print(het[, c("method", "Q", "Q_df", "Q_pval")], row.names = FALSE)
write.csv(het, "02.analysis/opengwas/mr_heterogeneity.csv", row.names = FALSE)
pleio <- mr_pleiotropy_test(dat)
print(pleio[, c("egger_intercept", "se", "pval")], row.names = FALSE)
write.csv(pleio, "02.analysis/opengwas/mr_pleiotropy.csv", row.names = FALSE)

# ---------- Steiger 方向性 ----------
cat("\n[5] Steiger 方向性检验\n")
dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
dat$r.outcome <- get_r_from_lor(dat$beta.outcome, dat$eaf.outcome,
                                ncase = 60801, ncontrol = 123504,   # CARDIoGRAMplusC4D
                                prevalence = 60801 / (60801 + 123504))
st <- tryCatch(directionality_test(dat), error = function(e) NULL)
if (!is.null(st)) {
  print(st[, c("correct_causal_direction", "steiger_pval")], row.names = FALSE)
  write.csv(st, "02.analysis/opengwas/mr_steiger.csv", row.names = FALSE)
}

# ---------- MR-PRESSO ----------
cat("\n[6] MR-PRESSO\n")
suppressMessages(library(MRPRESSO))
presso <- tryCatch(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                             SdOutcome = "se.outcome", SdExposure = "se.exposure",
                             OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
                             NbDistribution = 1000), error = function(e) paste("PRESSO 失败:", conditionMessage(e)))
if (is.character(presso)) cat("    ", presso, "\n") else {
  cat("    全局检验 P =", presso$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
  sink("02.analysis/opengwas/mr_presso.txt"); print(presso); sink()
  cat("    详情见 02.analysis/opengwas/mr_presso.txt\n")
}

# ---------- 绘图 ----------
cat("\n[7] 绘图\n")
p1 <- mr_scatter_plot(res, dat)
p2 <- mr_forest_plot(mr_singlesnp(dat))
p3 <- mr_funnel_plot(mr_singlesnp(dat))
p4 <- mr_leaveoneout_plot(mr_leaveoneout(dat))
ggplot2::ggsave("04.figures/opengwas_scatter.png", p1[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/opengwas_forest.png", p2[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/opengwas_funnel.png", p3[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/opengwas_leaveoneout.png", p4[[1]], width = 7, height = 8)
cat("    04.figures/opengwas_{scatter,forest,funnel,leaveoneout}.png 已保存\n")
cat("\nLDL-C -> CHD 分析完成 ✔\n")
