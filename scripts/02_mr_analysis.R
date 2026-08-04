#!/usr/bin/env Rscript
# 02_mr_analysis.R — MR 主分析：五种估计方法 + 异质性 + 多效性检验
# 输入: 02.analysis/harmonised.csv（由 01_data_prep.R 生成）
# 输出: 02.analysis/mr_results.csv、mr_heterogeneity.csv、mr_pleiotropy.csv
#       04.figures/mr_{scatter,forest,funnel}.png
# 方法: IVW / MR-Egger / Weighted median / Weighted mode / Simple mode

options(width = 150)
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("MR 主分析（Telomere_length -> CHD）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

dat <- read.csv("02.analysis/harmonised.csv", stringsAsFactors = FALSE)
dat <- subset(dat, mr_keep)
cat("[1] 工具变量数:", nrow(dat), "\n")

# ---------- 五种方法 ----------
cat("[2] 五种 MR 方法\n")
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                               "mr_weighted_median", "mr_weighted_mode",
                               "mr_simple_mode"))
res$OR <- exp(res$b)
res$OR_lci <- exp(res$b - 1.96 * res$se)
res$OR_uci <- exp(res$b + 1.96 * res$se)
print(res[, c("method", "nsnp", "b", "se", "pval", "OR", "OR_lci", "OR_uci")],
      row.names = FALSE)
write.csv(res, "02.analysis/mr_results.csv", row.names = FALSE)

# ---------- 异质性 & 多效性 ----------
cat("\n[3] 异质性（Cochran Q）与水平多效性（Egger 截距）\n")
het <- mr_heterogeneity(dat)
print(het[, c("method", "Q", "Q_df", "Q_pval")], row.names = FALSE)
write.csv(het, "02.analysis/mr_heterogeneity.csv", row.names = FALSE)

pleio <- mr_pleiotropy_test(dat)
print(pleio[, c("egger_intercept", "se", "pval")], row.names = FALSE)
write.csv(pleio, "02.analysis/mr_pleiotropy.csv", row.names = FALSE)

# ---------- 方向性检验（Steiger） ----------
cat("\n[4] Steiger 方向性检验\n")
# 暴露为连续性状（端粒长度），结局为二分类（CHD）
# 预计算 r 值: 暴露用 get_r_from_pn，结局用 get_r_from_lor（log OR）
dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
# CARDIoGRAM 2011 (Schunkert et al.): 22,233 病例 / 64,762 对照
dat$r.outcome <- get_r_from_lor(dat$beta.outcome, dat$eaf.outcome,
                                ncase = 22233, ncontrol = 64762,
                                prevalence = 22233 / (22233 + 64762))
steiger <- tryCatch(directionality_test(dat), error = function(e) {
  cat("    Steiger 检验失败:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(steiger)) {
  print(steiger[, c("snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction",
                    "steiger_pval")], row.names = FALSE)
  write.csv(steiger, "02.analysis/mr_steiger.csv", row.names = FALSE)
}

# ---------- 绘图 ----------
cat("\n[5] 绘图\n")
p1 <- mr_scatter_plot(res, dat)
p2 <- mr_forest_plot(mr_singlesnp(dat))
p3 <- mr_funnel_plot(mr_singlesnp(dat))
ggplot2::ggsave("04.figures/mr_scatter.png", p1[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/mr_forest.png", p2[[1]], width = 7, height = 6)
ggplot2::ggsave("04.figures/mr_funnel.png", p3[[1]], width = 7, height = 6)
cat("    已保存 04.figures/mr_{scatter,forest,funnel}.png\n")
cat("\n主分析完成 ✔\n")
