#!/usr/bin/env Rscript
# 48_drug_target_network.R — 药物-靶点网络分析
#
# 整合多蛋白靶点 pQTL/eQTL × 多结局，构建药物-靶点-疾病网络
#
# 靶点: PCSK9(pQTL) / HMGCR(eQTL) / CETP(eQTL) / NPC1L1(eQTL)
# 结局: CHD / CAD / MS / 白血病 / HOMA-IR / 雌二醇
# 输出: 网络图 + 矩阵表

options(width = 150)
dir.create("02.analysis/drug_network", showWarnings = FALSE)
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("药物-靶点网络分析\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 靶点与结局定义 ----------
# 靶点-结局文件映射（已有数据）
targets <- data.frame(
  id = c("PCSK9", "HMGCR", "CETP", "NPC1L1"),
  inst_file = c("02.analysis/opengwas/online_PCSK9_pqtl.csv",
                "02.analysis/multiomics/eqtl_HMGCR.csv",
                "02.analysis/drug_repurposing/inst_CETP.csv",
                "02.analysis/drug_repurposing/inst_NPC1L1.csv"),
  outcome_file = c("02.analysis/opengwas/online_outcome_PCSK9_CHD.csv",
                   "02.analysis/multiomics/outcome_HMGCR_CHD.csv",
                   "02.analysis/drug_repurposing/outcome_CETP_CHD.csv",
                   "02.analysis/drug_repurposing/outcome_NPC1L1_CHD.csv"),
  stringsAsFactors = FALSE)
outcomes <- data.frame(id = "ieu-a-7", trait = "CHD", stringsAsFactors = FALSE)

# ---------- IVW 批处理 ----------
mr_ivw_fast <- function(bx, by, sey) {
  if (length(bx) < 1) return(NA)
  w <- 1 / sey^2; b <- sum(bx * by * w) / sum(bx^2 * w)
  b
}

# 靶点-结局结果矩阵
res_matrix <- matrix(NA, nrow(targets), nrow(outcomes),
                     dimnames = list(targets$id, outcomes$trait))

for (t in seq_len(nrow(targets))) {
  cat(sprintf("\n[%s] 处理靶点: %s\n", targets$id[t], targets$id[t]))
  e <- tryCatch(read.table(targets$inst_file[t], header = TRUE, sep = "\t", stringsAsFactors = FALSE),
                error = function(e) data.frame())
  if (nrow(e) < 1) { cat("    无工具变量, 跳过\n"); next }
  cat("    工具变量:", nrow(e), "个\n")
  for (o in seq_len(nrow(outcomes))) {
    if (!file.exists(targets$outcome_file[t])) { cat("    ", outcomes$trait[o], ": 无结局文件\n"); next }
    o_data <- tryCatch(read.table(targets$outcome_file[t], header = TRUE, sep = "\t", stringsAsFactors = FALSE),
                      error = function(e) NULL)
    if (is.null(o_data) || nrow(o_data) < 1) { cat("    ", outcomes$trait[o], ": 空结局文件\n"); next }
    m <- merge(e, o_data, by = "rsid")
    if (nrow(m) < 1) { cat("    ", outcomes$trait[o], ": 无匹配 SNP\n"); next }
    m$flip <- m$ea.x != m$ea.y
    m$beta_y <- ifelse(m$flip, -m$beta.y, m$beta.y)
    b <- mr_ivw_fast(m$beta.x, m$beta_y, m$se.y)
    res_matrix[t, o] <- b
    cat(sprintf("    %s: beta=%.4f (n=%d)\n", outcomes$trait[o], b, nrow(m)))
  }
}

# ---------- 网络图 ----------
cat("\n\n绘图: 靶点-疾病网络热图\n")
res_long <- reshape2::melt(res_matrix, varnames = c("靶点", "结局"), value.name = "beta")
res_long <- res_long[!is.na(res_long$beta), ]
if (nrow(res_long) > 0) {
  res_long$效应方向 <- ifelse(res_long$beta > 0, "正向(危险)", "负向(保护)")
  p <- ggplot(res_long, aes(x = 结局, y = 靶点, fill = beta)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.3f", beta)), size = 3.5) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, name = "beta") +
    labs(title = "药物-靶点-疾病网络: 批量 MR 效应矩阵",
         subtitle = paste0(nrow(res_long), " 个靶点-结局对"),
         x = "", y = "药物靶点") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggplot2::ggsave("04.figures/drug_target_network.png", p, width = 8, height = 5)
  cat("    04.figures/drug_target_network.png 已保存\n")
  write.csv(res_long, "02.analysis/drug_network/network_matrix.csv", row.names = FALSE)
  print(res_long, row.names = FALSE)
} else {
  cat("    无结果可绘图\n")
}
cat("\n药物-靶点网络分析完成 ✔\n")