#!/usr/bin/env Rscript
# 24_opengwas_multi.R — OpenGWAS 在线数据: LDL/HDL/TG -> CHD 多暴露对比
#
# 输入（由 22_opengwas_api.py 在线获取）:
#   暴露: online_instruments_ieu-b-1{10,09,11}.csv (tophits + 本地 clumping)
#   结局: online_outcome_ieu-a-7 对应查询结果
# 输出: 02.analysis/opengwas/multi/ 对比表 + 04.figures/ 森林图

options(width = 150)
dir.create("02.analysis/opengwas/multi", showWarnings = FALSE)
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("OpenGWAS 在线数据多暴露对比: LDL/HDL/TG -> CHD\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

exposures <- data.frame(
  id = c("ieu-b-110", "ieu-b-109", "ieu-b-111"),
  trait = c("LDL cholesterol", "HDL cholesterol", "Triglycerides"),
  n = c(188578, 188578, 188578), stringsAsFactors = FALSE)

all_res <- list()
for (i in seq_len(nrow(exposures))) {
  ex <- exposures[i, ]
  cat("----------------------------------------\n")
  cat(sprintf("[%d] %s (%s)\n", i, ex$trait, ex$id))
  expo_raw <- read.table(sprintf("02.analysis/opengwas/online_instruments_%s.csv", ex$id),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  outc_raw <- read.table(sprintf("02.analysis/opengwas/online_outcome_%s.csv", ex$id),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  expo_f <- format_data(expo_raw, type = "exposure", header = TRUE,
                        snp_col = "rsid", beta_col = "beta", se_col = "se",
                        effect_allele_col = "ea", other_allele_col = "nea",
                        eaf_col = "eaf", pval_col = "p")
  outc_f <- format_data(outc_raw, type = "outcome", header = TRUE,
                        snp_col = "rsid", beta_col = "beta", se_col = "se",
                        effect_allele_col = "ea", other_allele_col = "nea",
                        eaf_col = "eaf", pval_col = "p")
  dat <- subset(harmonise_data(expo_f, outc_f), mr_keep)
  cat("    harmonise 后工具变量:", nrow(dat), "\n")

  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                                 "mr_weighted_mode", "mr_simple_mode"))
  res$trait <- ex$trait
  res$exposure_id <- ex$id
  res$OR <- exp(res$b); res$OR_lci <- exp(res$b - 1.96 * res$se); res$OR_uci <- exp(res$b + 1.96 * res$se)
  print(res[, c("method", "nsnp", "b", "se", "pval", "OR", "OR_lci", "OR_uci")], row.names = FALSE)

  het <- mr_heterogeneity(dat); pleio <- mr_pleiotropy_test(dat)
  res$Q_pval <- ifelse(nrow(het) > 0, het$Q_pval[1], NA)
  res$egger_intercept_p <- ifelse(nrow(pleio) > 0, pleio$pval[1], NA)
  all_res[[i]] <- res
  write.csv(dat, sprintf("02.analysis/opengwas/multi/harmonised_%s.csv", ex$id), row.names = FALSE)
}

# ---------- 汇总对比表 ----------
cat("\n========================================\n")
cat("多暴露对比汇总（IVW 与加权中位数）\n")
res_all <- do.call(rbind, all_res)
summ <- subset(res_all, method %in% c("Inverse variance weighted", "Weighted median"))
summ <- summ[, c("trait", "method", "nsnp", "b", "se", "pval", "OR", "OR_lci", "OR_uci",
                 "Q_pval", "egger_intercept_p")]
print(summ, row.names = FALSE)
write.csv(summ, "02.analysis/opengwas/multi/comparison_summary.csv", row.names = FALSE)

# ---------- 森林图 ----------
cat("\n绘图: 多暴露森林图\n")
summ$label <- paste0(summ$trait, " (", summ$method, ")")
p <- ggplot(summ, aes(x = OR, y = reorder(label, OR), color = trait)) +
  geom_point(size = 3) + geom_errorbarh(aes(xmin = OR_lci, xmax = OR_uci), height = 0.2) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10() +
  labs(title = "LDL/HDL/TG -> CHD 多暴露 MR 对比 (OR 与 95%CI)", x = "OR (log scale)", y = "") +
  theme_bw()
ggplot2::ggsave("04.figures/online_multi_forest.png", p, width = 8, height = 5)
cat("    04.figures/online_multi_forest.png 已保存\n")
cat("\n多暴露对比完成 ✔\n")
