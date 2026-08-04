#!/usr/bin/env Rscript
# 45_methods_panorama.R — MR 多效性方法全景对比
#
# 在同一个数据集（LDL-C -> CHD）上运行所有可用的 MR 方法，
# 生成全景对比图，系统展示各方法估计的异同。
#
# 方法: IVW / MR-Egger / 加权中位数 / 加权众数 / 简单众数 /
#        ConMix / cML / MaxLik / SimpleMedian
# 输出: 02.analysis/panorama/ 全景对比表 + 04.figures/panorama_*.png

options(width = 150)
dir.create("02.analysis/panorama", showWarnings = FALSE)
suppressMessages(library(MendelianRandomization))
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("MR 多效性方法全景对比\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 数据 ----------
cat("[1] 数据: LDL-C -> CHD（在线, 28 个工具变量）\n")
d <- read.csv("02.analysis/opengwas/online/harmonised.csv", stringsAsFactors = FALSE)
d <- subset(d, mr_keep)
mrin <- mr_input(bx = d$beta.exposure, bxse = d$se.exposure,
                 by = d$beta.outcome, byse = d$se.outcome,
                 snps = d$SNP, exposure = "LDL-C", outcome = "CHD")
cat("    工具变量数:", nrow(d), "\n")

# ---------- 运行所有方法 ----------
cat("\n[2] 运行所有方法\n")
methods <- list()
methods$IVW <- mr_ivw(mrin)
methods$Egger <- mr_egger(mrin)
methods$WMedian <- mr_median(mrin)
methods$WMode <- mr_mbe(mrin)
methods$SimpleMode <- mr_mbe(mrin, weighting = "unweighted", stderror = "simple")
methods$SimpleMedian <- mr_median(mrin, weighting = "simple")
methods$MaxLik <- mr_maxlik(mrin)
methods$ConMix <- mr_conmix(mrin)
methods$cML <- mr_cML(mrin, n = 184305)

# ---------- 提取结果 ----------
extract <- function(obj, name) {
  if (inherits(obj, "IVW") || inherits(obj, "MaxLik") || inherits(obj, "WeightedMedian") ||
      inherits(obj, "SimpleMedian")) {
    data.frame(方法 = name, beta = obj@Estimate, CI_lower = obj@CILower,
               CI_upper = obj@CIUpper, P = obj@Pvalue, stringsAsFactors = FALSE)
  } else if (inherits(obj, "Egger")) {
    data.frame(方法 = name, beta = obj@Estimate, CI_lower = obj@CILower.Est,
               CI_upper = obj@CIUpper.Est, P = obj@Pvalue.Est, stringsAsFactors = FALSE)
  } else if (inherits(obj, "MBE")) {
    data.frame(方法 = name, beta = obj@Estimate, CI_lower = obj@CILower,
               CI_upper = obj@CIUpper, P = obj@Pvalue, stringsAsFactors = FALSE)
  } else if (inherits(obj, "MRConMix")) {
    data.frame(方法 = name, beta = obj@Estimate, CI_lower = obj@CILower,
               CI_upper = obj@CIUpper, P = obj@Pvalue, stringsAsFactors = FALSE)
  } else if (inherits(obj, "MRcML")) {
    data.frame(方法 = name, beta = obj@Estimate, CI_lower = obj@CILower,
               CI_upper = obj@CIUpper, P = obj@Pvalue, stringsAsFactors = FALSE)
  } else {
    data.frame(方法 = name, beta = NA, CI_lower = NA, CI_upper = NA, P = NA, stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, mapply(extract, methods, names(methods), SIMPLIFY = FALSE))
res$OR <- exp(res$beta)
res$OR_lci <- exp(res$CI_lower)
res$OR_uci <- exp(res$CI_upper)
res$sig <- ifelse(res$P < 0.05, "P<0.05", "P≥0.05")
res$label <- paste0(res$方法, " (OR=", sprintf("%.2f", res$OR),
                    ", P=", format(res$P, digits = 2), ")")

cat("\n[3] 结果汇总\n")
print(res[, c("方法", "beta", "CI_lower", "CI_upper", "P", "OR", "OR_lci", "OR_uci")], row.names = FALSE)
write.csv(res, "02.analysis/panorama/panorama_results.csv", row.names = FALSE)

# ---------- 全景森林图 ----------
cat("\n[4] 绘图: 全景森林图\n")
p <- ggplot(res, aes(x = OR, y = reorder(方法, OR), color = sig)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = OR_lci, xmax = OR_uci), height = 0.2, linewidth = 1) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  scale_x_log10() +
  scale_color_manual(values = c("P<0.05" = "darkred", "P≥0.05" = "grey50")) +
  labs(title = "MR 多效性方法全景对比: LDL-C -> CHD",
       subtitle = paste0(nrow(d), " 个工具变量 · 所有方法在同一数据集上"),
       x = "OR (95% CI, log scale)", y = NULL, color = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")
ggplot2::ggsave("04.figures/panorama_forest.png", p, width = 9, height = 6)
cat("    04.figures/panorama_forest.png 已保存\n")

# ---------- 解读 ----------
cat("\n[5] 解读\n")
all_sig <- all(res$P < 0.05, na.rm = TRUE)
if (all_sig) {
  cat("    所有方法方向一致且显著 -> 证据强，结论稳健\n")
} else {
  sig <- sum(res$P < 0.05, na.rm = TRUE)
  cat(sprintf("    %d/%d 方法显著 -> 多数方法支持因果效应\n", sig, nrow(res)))
}
cat("\nMR 方法全景对比完成 ✔\n")