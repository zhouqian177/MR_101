#!/usr/bin/env Rscript
# 38_drug_repurposing.R — 药物重定位扫描（Drug Repurposing MR）
#
# 思路: 用多个已知药物靶点的 pQTL/eQTL 作为工具变量, 对疾病结局批量 MR,
#       识别"调节该靶点可降低疾病风险"的证据 -> 潜在可重定位的药物靶点。
#
# 靶点（均为降脂/心血管药物靶点）:
#   - PCSK9  (蛋白, pQTL)   -> PCSK9 抑制剂(依洛尤单抗等)
#   - HMGCR  (基因表达)     -> 他汀类
#   - CETP   (基因表达)     -> CETP 抑制剂(安塞曲匹等)
#   - NPC1L1 (基因表达)     -> 依折麦布
#
# 结局: CHD (ieu-a-7)
# 输出: 02.analysis/drug_repurposing/ 扫描结果 + 04.figures/repurposing_forest.png

options(width = 150)
dir.create("02.analysis/drug_repurposing", showWarnings = FALSE)
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("药物重定位扫描（多靶点 -> CHD）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 通用: 暴露/结局合并 + IVW ----------
mr_pair <- function(expo_f, outc_f, target, drug, omics) {
  e <- read.table(expo_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  o <- read.table(outc_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(e) == 0 || nrow(o) == 0) return(NULL)
  m <- merge(e, o, by = "rsid")
  if (nrow(m) == 0) return(NULL)
  m$flip <- m$ea.x != m$ea.y
  m$beta_y <- ifelse(m$flip, -m$beta.y, m$beta.y)
  m$beta_x <- m$beta.x; m$se_y <- m$se.y
  w <- 1 / m$se_y^2
  b <- sum(m$beta_x * m$beta_y * w) / sum(m$beta_x^2 * w)
  se <- sqrt(1 / sum(m$beta_x^2 * w))
  p <- 2 * pnorm(-abs(b / se))
  data.frame(靶点 = target, 药物 = drug, 组学 = omics, 工具数 = nrow(m),
             beta = b, SE = se, P = p, OR = exp(b),
             OR_lci = exp(b - 1.96 * se), OR_uci = exp(b + 1.96 * se))
}

# ---------- 1. 各靶点 MR ----------
cat("[1] PCSK9（蛋白, pQTL）-> CHD\n")
r1 <- mr_pair("02.analysis/opengwas/online_PCSK9_pqtl.csv",
              "02.analysis/opengwas/online_outcome_PCSK9_CHD.csv",
              "PCSK9", "PCSK9 抑制剂(单抗/ASO)", "pQTL")
if (!is.null(r1)) print(r1, row.names = FALSE)

cat("\n[2] HMGCR（基因表达, eQTL）-> CHD\n")
r2 <- mr_pair("02.analysis/multiomics/eqtl_HMGCR.csv",
              "02.analysis/multiomics/outcome_HMGCR_CHD.csv",
              "HMGCR", "他汀类", "eQTL")
if (!is.null(r2)) print(r2, row.names = FALSE)

cat("\n[3] CETP（基因表达, eQTL）-> CHD\n")
r3 <- mr_pair("02.analysis/drug_repurposing/inst_CETP.csv",
              "02.analysis/drug_repurposing/outcome_CETP_CHD.csv",
              "CETP", "CETP 抑制剂(他塞曲匹)", "eQTL")
if (!is.null(r3)) print(r3, row.names = FALSE)

cat("\n[4] NPC1L1（基因表达, eQTL）-> CHD\n")
r4 <- mr_pair("02.analysis/drug_repurposing/inst_NPC1L1.csv",
              "02.analysis/drug_repurposing/outcome_NPC1L1_CHD.csv",
              "NPC1L1", "依折麦布", "eQTL")
if (!is.null(r4)) print(r4, row.names = FALSE)

# ---------- 2. 扫描汇总 ----------
cat("\n[5] 药物重定位扫描汇总\n")
res <- do.call(rbind, Filter(Negate(is.null), list(r1, r2, r3, r4)))
res$sig <- ifelse(res$P < 0.05, "P<0.05", "NS")
print(res, row.names = FALSE)
write.csv(res, "02.analysis/drug_repurposing/repurposing_scan.csv", row.names = FALSE)

# ---------- 3. 森林图 ----------
cat("\n[6] 绘图: 药物靶点重定位森林图\n")
res$label <- paste0(res$靶点, " (", res$组学, ", n=", res$工具数, ")")
p <- ggplot(res, aes(x = OR, y = reorder(label, OR), color = sig)) +
  geom_point(size = 4) + geom_errorbarh(aes(xmin = OR_lci, xmax = OR_uci), height = 0.2) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10() +
  scale_color_manual(values = c("P<0.05" = "darkred", "NS" = "grey50")) +
  labs(title = "药物重定位扫描: 药物靶点 -> CHD (IVW)", x = "OR (log)", y = "") +
  theme_bw()
ggplot2::ggsave("04.figures/repurposing_forest.png", p, width = 7, height = 4)
cat("    04.figures/repurposing_forest.png 已保存\n")

# ---------- 4. 解读 ----------
cat("\n[7] 解读\n")
sig_hits <- res[res$P < 0.05, ]
if (nrow(sig_hits) > 0) {
  cat("    显著靶点（P<0.05）: ", paste(sig_hits$靶点, collapse = ", "), "\n")
  for (i in seq_len(nrow(sig_hits))) {
    cat(sprintf("      %s: OR=%.3f (%.3f-%.3f) -> 提示 %s 具有降低 CHD 风险潜力（重定位候选）\n",
                sig_hits$靶点[i], sig_hits$OR[i], sig_hits$OR_lci[i], sig_hits$OR_uci[i],
                sig_hits$药物[i]))
  }
} else {
  cat("    无显著靶点; 或工具变量功效不足（单/少数 cis 变异）\n")
}
cat("\n药物重定位扫描完成 ✔\n")
