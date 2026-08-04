#!/usr/bin/env Rscript
# 34_mr_scan.R — 全基因组 MR 扫描（MR Scan）：LDL-C 工具变量 × 多结局批量 MR
#
# 数据（OpenGWAS 在线）:
#   暴露: LDL-C (ieu-b-110), 28 个工具变量（已 clumping）
#   结局: 6 个数据集（冠心病×2、CAD、多发性硬化、淋巴白血病、HOMA-IR、雌二醇）
# 输出: 02.analysis/mr_scan/ 扫描结果表 + 04.figures/mr_scan_*.png

options(width = 150)
dir.create("02.analysis/mr_scan", showWarnings = FALSE)
suppressMessages(library(TwoSampleMR))
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("全基因组 MR 扫描（LDL-C × 6 结局）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

outcomes <- data.frame(
  id = c("ieu-a-7", "ieu-a-8", "ebi-a-GCST005194", "ieu-b-18", "ieu-b-4956", "ieu-b-118", "ieu-b-4872"),
  trait = c("CHD (CARDIoGRAM+C4D)", "CHD (CARDIoGRAM 2011)", "CAD (2022)",
            "Multiple sclerosis", "Lymphoid leukaemia", "HOMA-IR (T2D 相关)", "Oestradiol"),
  stringsAsFactors = FALSE)

# 暴露数据（在线 harmonised 数据中的 28 个工具变量）
d <- read.csv("02.analysis/opengwas/online/harmonised.csv", stringsAsFactors = FALSE)
expo <- d[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
              "beta.exposure", "se.exposure", "pval.exposure", "eaf.exposure")]
colnames(expo) <- c("SNP", "ea", "nea", "beta", "se", "p", "eaf")

scan <- list()
for (i in seq_len(nrow(outcomes))) {
  oc <- outcomes[i, ]
  f <- sprintf("02.analysis/mr_scan/outcome_%s.csv", oc$id)
  if (!file.exists(f)) { cat("跳过缺失:", oc$id, "\n"); next }
  outc <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat(sprintf("[%d] %s (%s): %d 条关联\n", i, oc$trait, oc$id, nrow(outc)))
  if (nrow(outc) < 3) { cat("    SNP 过少, 跳过\n"); next }

  expo_f <- format_data(expo, type = "exposure", header = TRUE,
                        snp_col = "SNP", beta_col = "beta", se_col = "se",
                        effect_allele_col = "ea", other_allele_col = "nea",
                        eaf_col = "eaf", pval_col = "p")
  outc_f <- format_data(outc, type = "outcome", header = TRUE,
                        snp_col = "rsid", beta_col = "beta", se_col = "se",
                        effect_allele_col = "ea", other_allele_col = "nea",
                        eaf_col = "eaf", pval_col = "p")
  dat <- subset(harmonise_data(expo_f, outc_f), mr_keep)
  if (nrow(dat) < 3) { cat("    harmonise 后 SNP 过少, 跳过\n"); next }

  res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median"))
  res$trait <- oc$trait; res$outcome_id <- oc$id
  res$OR <- exp(res$b); res$OR_lci <- exp(res$b - 1.96 * res$se); res$OR_uci <- exp(res$b + 1.96 * res$se)
  scan[[oc$id]] <- res
  cat(sprintf("    IVW OR=%.2f (%.2f-%.2f, P=%.2g) | WMedian OR=%.2f (P=%.2g)\n",
              res$OR[1], res$OR_lci[1], res$OR_uci[1], res$pval[1],
              exp(res$b[2]), res$pval[2]))
}

scan_all <- do.call(rbind, scan)
scan_all$sig <- ifelse(scan_all$pval < 0.05, "P<0.05", "NS")
write.csv(scan_all, "02.analysis/mr_scan/scan_results.csv", row.names = FALSE)
cat("\n========================================\n")
cat("MR 扫描汇总（IVW）\n")
ivw_res <- subset(scan_all, method == "Inverse variance weighted")
print(ivw_res[, c("trait", "nsnp", "b", "se", "pval", "OR", "OR_lci", "OR_uci")], row.names = FALSE)

# ---------- 森林图（批量） ----------
cat("\n绘图: MR 扫描森林图\n")
ivw_res$label <- paste0(ivw_res$trait, " (n=", ivw_res$nsnp, ")")
p <- ggplot(ivw_res, aes(x = OR, y = reorder(label, OR), color = sig)) +
  geom_point(size = 3) + geom_errorbarh(aes(xmin = OR_lci, xmax = OR_uci), height = 0.2) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10() +
  scale_color_manual(values = c("P<0.05" = "darkred", "NS" = "grey50")) +
  labs(title = "MR 扫描: LDL-C 对多结局的因果效应 (IVW)", x = "OR (log)", y = "") +
  theme_bw()
ggplot2::ggsave("04.figures/mr_scan_forest.png", p, width = 8, height = 5)
cat("    04.figures/mr_scan_forest.png 已保存\n")
cat("\nMR 扫描完成 ✔\n")
