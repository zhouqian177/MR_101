#!/usr/bin/env Rscript
# 49_gwas_mr_scan.R — 全基因组 MR 批量扫描（多暴露×多结局矩阵）
#
# 扩展 34_mr_scan.R 为多暴露×多结局矩阵：
# - 暴露: LDL-C / HDL-C / TG（已有数据，OpenGWAS 在线）
# - 结局: CHD / CAD / MS / 白血病 / HOMA-IR / 雌二醇
# 输出: 效应矩阵表 + 热图

options(width = 150)
dir.create("02.analysis/gwas_scan", showWarnings = FALSE)
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

cat("========================================\n")
cat("全基因组 MR 批量扫描\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 暴露定义 ----------
exposures <- data.frame(
  id = c("ieu-b-110", "ieu-b-109", "ieu-b-111"),
  trait = c("LDL-C", "HDL-C", "TG"),
  inst_file = c("02.analysis/opengwas/online_instruments_ieu-b-110.csv",
                "02.analysis/opengwas/online_instruments_ieu-b-109.csv",
                "02.analysis/opengwas/online_instruments_ieu-b-111.csv"),
  outcome_file = c("02.analysis/opengwas/online_outcome_ieu-b-110.csv",
                   "02.analysis/opengwas/online_outcome_ieu-b-109.csv",
                   "02.analysis/opengwas/online_outcome_ieu-b-111.csv"),
  stringsAsFactors = FALSE)

# ---------- IVW 函数 ----------
ivw_fast <- function(bx, by, sey) {
  if (length(bx) < 3) return(NA)
  w <- 1 / sey^2; b <- sum(bx * by * w) / sum(bx^2 * w)
  b
}

# ---------- 批量 MR ----------
res <- data.frame()
for (e_idx in seq_len(nrow(exposures))) {
  ex <- exposures[e_idx, ]
  cat(sprintf("\n[%d/%d] 暴露: %s (%s)\n", e_idx, nrow(exposures), ex$trait, ex$id))
  e <- tryCatch(read.table(ex$inst_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
                error = function(e) data.frame())
  if (nrow(e) < 1) { cat("    无工具变量, 跳过\n"); next }
  # 对每个结局文件做 MR
  o_files <- list.files("02.analysis/mr_scan", pattern = "outcome_.*\\.csv", full.names = TRUE)
  # 也包含每个暴露对应的 outcome 文件
  if (file.exists(ex$outcome_file)) o_files <- c(o_files, ex$outcome_file)
  o_files <- unique(o_files)
  for (of in o_files) {
    o_data <- tryCatch(read.table(of, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
                       error = function(e) NULL)
    if (is.null(o_data) || nrow(o_data) < 1) next
    # 识别结局名称
    oid <- gsub(".*outcome_|_CHD\\.csv|_ieu-a-7\\.csv|_ieu-b-.*\\.csv", "", basename(of))
    oid <- gsub("ieu-a-7|ieu-b-[0-9]+", "", oid)
    oid <- ifelse(nchar(oid) < 2, gsub("outcome_|_ieu-.*\\.csv", "", basename(of)), oid)
    m <- merge(e, o_data, by = "rsid")
    if (nrow(m) < 3) next
    m$flip <- if ("ea.x" %in% colnames(m) && "ea.y" %in% colnames(m)) m$ea.x != m$ea.y else FALSE
    m$beta_y <- ifelse(m$flip, -m$beta.y, m$beta.y)
    b <- ivw_fast(m$beta.x, m$beta_y, m$se.y)
    if (!is.na(b)) {
      res <- rbind(res, data.frame(暴露 = ex$trait, 结局文件 = basename(of), beta = b,
                                    n = nrow(m), stringsAsFactors = FALSE))
    }
  }
}
cat(sprintf("\n[完成] 共 %d 个暴露-结局对\n", nrow(res)))

# ---------- 热图 ----------
if (nrow(res) > 0) {
  res$结局标签 <- gsub("\\.csv", "", gsub("outcome_", "", res$结局文件))
  res$效应方向 <- ifelse(res$beta > 0, "正向", "负向")
  p <- ggplot(res, aes(x = 结局标签, y = 暴露, fill = beta)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.3f", beta)), size = 3.5) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, name = "beta") +
    labs(title = "全基因组 MR 批量扫描: 多暴露×多结局",
         x = "结局", y = "暴露") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggplot2::ggsave("04.figures/gwas_scan_heatmap.png", p, width = 8, height = 5)
  cat("    04.figures/gwas_scan_heatmap.png 已保存\n")
  write.csv(res, "02.analysis/gwas_scan/scan_matrix.csv", row.names = FALSE)
  print(res, row.names = FALSE)
}
cat("\n全基因组 MR 批量扫描完成 ✔\n")