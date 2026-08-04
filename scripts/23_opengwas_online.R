#!/usr/bin/env Rscript
# 23_opengwas_online.R — OpenGWAS 在线数据: LDL-C -> CHD 全面 MR 分析
#
# 输入（由 22_opengwas_api.py 从 OpenGWAS API 在线获取）:
#   暴露: 02.analysis/opengwas/online_instruments_ldl.csv (tophits P<5e-8 + 本地 clumping)
#   结局: 02.analysis/opengwas/online_outcome_chd.csv (associations 查询)
# 输出: 02.analysis/opengwas/online_*/ 结果表与 04.figures/online_*.png
#
# 说明: API 数据获取用 Python 通道(R libcurl 过旧), 分析层用 R(TwoSampleMR)

options(width = 150)
dir.create("02.analysis/opengwas/online", showWarnings = FALSE)
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("OpenGWAS 在线数据 MR: LDL-C -> CHD\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 读取在线数据 ----------
expo_raw <- read.table("02.analysis/opengwas/online_instruments_ldl.csv",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
outc_raw <- read.table("02.analysis/opengwas/online_outcome_chd.csv",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("[1] 在线数据: 暴露(ieu-b-110 LDL) 工具变量", nrow(expo_raw),
    "| 结局(ieu-a-7 CHD) 关联", nrow(outc_raw), "\n")

# ---------- 2. harmonise ----------
cat("[2] harmonise（统一效应等位基因）\n")
expo_f <- format_data(expo_raw, type = "exposure", header = TRUE,
                      snp_col = "rsid", beta_col = "beta", se_col = "se",
                      effect_allele_col = "ea", other_allele_col = "nea",
                      eaf_col = "eaf", pval_col = "p")
outc_f <- format_data(outc_raw, type = "outcome", header = TRUE,
                      snp_col = "rsid", beta_col = "beta", se_col = "se",
                      effect_allele_col = "ea", other_allele_col = "nea",
                      eaf_col = "eaf", pval_col = "p")
dat <- harmonise_data(expo_f, outc_f)
dat <- subset(dat, mr_keep)
cat("    harmonise 后工具变量:", nrow(dat), "\n")
write.csv(dat, "02.analysis/opengwas/online/harmonised.csv", row.names = FALSE)

# ---------- 3. 五种方法 ----------
cat("\n[3] 五种 MR 方法\n")
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                               "mr_weighted_mode", "mr_simple_mode"))
res$OR <- exp(res$b); res$OR_lci <- exp(res$b - 1.96 * res$se); res$OR_uci <- exp(res$b + 1.96 * res$se)
print(res[, c("method", "nsnp", "b", "se", "pval", "OR", "OR_lci", "OR_uci")], row.names = FALSE)
write.csv(res, "02.analysis/opengwas/online/mr_results.csv", row.names = FALSE)

# ---------- 4. 异质性 / 多效性 ----------
cat("\n[4] 异质性与水平多效性\n")
het <- mr_heterogeneity(dat)
print(het[, c("method", "Q", "Q_df", "Q_pval")], row.names = FALSE)
write.csv(het, "02.analysis/opengwas/online/mr_heterogeneity.csv", row.names = FALSE)
pleio <- mr_pleiotropy_test(dat)
print(pleio[, c("egger_intercept", "se", "pval")], row.names = FALSE)
write.csv(pleio, "02.analysis/opengwas/online/mr_pleiotropy.csv", row.names = FALSE)

# ---------- 5. Steiger 方向性 ----------
cat("\n[5] Steiger 方向性检验\n")
dat$r.exposure <- get_r_from_pn(dat$pval.exposure, 188578)
dat$r.outcome <- get_r_from_lor(dat$beta.outcome, dat$eaf.outcome,
                                ncase = 60801, ncontrol = 123504,
                                prevalence = 60801 / (60801 + 123504))
st <- tryCatch(directionality_test(dat), error = function(e) NULL)
if (!is.null(st)) {
  print(st[, c("correct_causal_direction", "steiger_pval")], row.names = FALSE)
  write.csv(st, "02.analysis/opengwas/online/mr_steiger.csv", row.names = FALSE)
}

# ---------- 6. MR-PRESSO ----------
cat("\n[6] MR-PRESSO\n")
suppressMessages(library(MRPRESSO))
presso <- tryCatch(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                             SdOutcome = "se.outcome", SdExposure = "se.exposure",
                             OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat,
                             NbDistribution = 1000), error = function(e) paste("PRESSO 失败:", conditionMessage(e)))
if (is.character(presso)) cat("    ", presso, "\n") else {
  cat("    全局检验 P =", presso$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
  sink("02.analysis/opengwas/online/mr_presso.txt"); print(presso); sink()
  cat("    详情见 02.analysis/opengwas/online/mr_presso.txt\n")
}

# ---------- 7. 绘图 ----------
cat("\n[7] 绘图\n")
p1 <- mr_scatter_plot(res, dat)
p2 <- mr_forest_plot(mr_singlesnp(dat))
p3 <- mr_funnel_plot(mr_singlesnp(dat))
p4 <- mr_leaveoneout_plot(mr_leaveoneout(dat))
ggplot2::ggsave("04.figures/online_scatter.png", p1[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/online_forest.png", p2[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/online_funnel.png", p3[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/online_leaveoneout.png", p4[[1]], width = 7, height = 8)
cat("    04.figures/online_{scatter,forest,funnel,leaveoneout}.png 已保存\n")
cat("\nOpenGWAS 在线 LDL-C -> CHD 分析完成 ✔\n")
