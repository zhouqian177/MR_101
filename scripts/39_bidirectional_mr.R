#!/usr/bin/env Rscript
# 39_bidirectional_mr.R — 双向 MR（Bidirectional MR）
#
# 背景: 标准 MR 假设暴露->结局因果方向; 但有些关联是双向的或反向因果。
#       双向 MR 在两个方向上都做两样本 MR, 对比效应大小与显著性。
#       正向显著 + 反向不显著 = 支持单向因果。
#       双向都显著 = 双向因果或反馈环路。
#
# 实例: LDL-C (ieu-b-110) ↔ CHD (ieu-a-7)
#   正向: LDL-C -> CHD（已有 results/mr_results.csv）
#   反向: CHD  -> LDL-C（本脚本）
# 输出: 02.analysis/bidirectional/ 对比表 + 04.figures/bidirectional_*.png

options(width = 150)
dir.create("02.analysis/bidirectional", showWarnings = FALSE)
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("双向 MR: LDL-C ↔ CHD\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 正向（已有结果）----------
cat("[1] 正向 MR: LDL-C -> CHD\n")
fwd <- read.csv("02.analysis/opengwas/online/mr_results.csv", stringsAsFactors = FALSE)
f_ivw <- subset(fwd, method == "Inverse variance weighted")
cat(sprintf("    IVW: beta=%.4f (SE=%.4f, P=%.3g), OR=%.3f\n",
            f_ivw$b, f_ivw$se, f_ivw$pval, exp(f_ivw$b)))

# ---------- 反向 MR: CHD -> LDL-C ----------
cat("\n[2] 反向 MR: CHD -> LDL-C\n")
e <- read.table("02.analysis/bidirectional/inst_CHD.csv", header = TRUE,
                sep = "\t", stringsAsFactors = FALSE)
o <- read.table("02.analysis/bidirectional/outcome_CHD_LDL.csv", header = TRUE,
                sep = "\t", stringsAsFactors = FALSE)
m <- merge(e, o, by = "rsid")
m$flip <- m$ea.x != m$ea.y
m$beta_y <- ifelse(m$flip, -m$beta.y, m$beta.y)
m$beta_x <- m$beta.x; m$se_y <- m$se.y
w <- 1 / m$se_y^2
b_rev <- sum(m$beta_x * m$beta_y * w) / sum(m$beta_x^2 * w)
se_rev <- sqrt(1 / sum(m$beta_x^2 * w))
p_rev <- 2 * pnorm(-abs(b_rev / se_rev))
cat(sprintf("    IVW: beta=%.4f (SE=%.4f, P=%.3g)\n", b_rev, se_rev, p_rev))
cat(sprintf("    工具变量数: %d\n", nrow(m)))

# ---------- 双向对比 ----------
cat("\n[3] 双向对比\n")
comp <- data.frame(
  方向 = c("LDL-C -> CHD", "CHD -> LDL-C"),
  beta = c(f_ivw$b, b_rev),
  SE = c(f_ivw$se, se_rev),
  P = c(f_ivw$pval, p_rev),
  OR = c(exp(f_ivw$b), exp(b_rev)),
  OR_lci = c(exp(f_ivw$b - 1.96 * f_ivw$se), exp(b_rev - 1.96 * se_rev)),
  OR_uci = c(exp(f_ivw$b + 1.96 * f_ivw$se), exp(b_rev + 1.96 * se_rev)))
print(comp, row.names = FALSE)
write.csv(comp, "02.analysis/bidirectional/comparison.csv", row.names = FALSE)

# 解读
if (p_rev >= 0.05) {
  cat("\n结论: 反向 (CHD->LDL) 不显著 -> 支持 LDL-C 对 CHD 的单向因果效应\n")
} else {
  cat("\n结论: 双向皆显著 -> 可能存在双向因果或反馈环路\n")
}

# ---------- 图 ----------
cat("\n[4] 绘图: 双向森林图\n")
comp$label <- paste0(comp$方向, " (P=", format(comp$P, digits = 3), ")")
p <- ggplot(comp, aes(x = OR, y = reorder(label, OR))) +
  geom_point(size = 4, color = "darkred") +
  geom_errorbarh(aes(xmin = OR_lci, xmax = OR_uci), height = 0.2, color = "darkred") +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10() +
  labs(title = "双向 MR: LDL-C ↔ CHD (IVW)", x = "OR (log)", y = "") +
  theme_bw()
ggplot2::ggsave("04.figures/bidirectional_forest.png", p, width = 6, height = 3)
cat("    04.figures/bidirectional_forest.png 已保存\n")
cat("\n双向 MR 完成 ✔\n")