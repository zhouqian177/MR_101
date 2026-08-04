#!/usr/bin/env Rscript
# 35_visual_report.R — MR 结果可视化报告（批量森林图 + 热图）
#
# 汇总本仓库各模块分析结果, 生成统一可视化:
#   1) 多模块森林图: 各分析方法 OR(95%CI) 汇总
#   2) 多暴露-多结局热图: LDL/HDL/TG x 结局的因果效应矩阵
#   3) 可视化报告总览（表格）
#
# 输入: 02.analysis/ 下各模块结果 CSV
# 输出: 02.analysis/visual/ 表格 + 04.figures/visual_*.png

options(width = 150)
dir.create("02.analysis/visual", showWarnings = FALSE)
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

cat("========================================\n")
cat("MR 结果可视化报告\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 收集各模块 MR 结果 ----------
cat("[1] 汇总各模块结果\n")
collect <- function(f, trait_label, source) {
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (!"method" %in% colnames(d)) return(NULL)
  d$trait <- trait_label
  d$source <- source
  d
}
parts <- list(
  collect("02.analysis/mr_results.csv", "Telomere_length->CHD", "two-sample"),
  collect("02.analysis/opengwas/online/mr_results.csv", "LDL-C->CHD(online)", "opengwas"),
  collect("02.analysis/mr_scan/scan_results.csv", "MR scan", "scan"))
parts <- parts[!sapply(parts, is.null)]
# 各模块结果列不同, 只保留共同列后再合并
common_cols <- Reduce(intersect, lapply(parts, colnames))
res_all <- do.call(rbind, lapply(parts, function(d) d[, common_cols]))
res_all$OR <- ifelse(is.na(res_all$OR), exp(res_all$b), res_all$OR)
res_all$OR_lci <- ifelse(is.na(res_all$OR_lci), exp(res_all$b - 1.96 * res_all$se), res_all$OR_lci)
res_all$OR_uci <- ifelse(is.na(res_all$OR_uci), exp(res_all$b + 1.96 * res_all$se), res_all$OR_uci)
res_all$sig <- ifelse(res_all$pval < 0.05, "P<0.05", "NS")
cat("    汇总结果行数:", nrow(res_all), "\n")

# ---------- 2. 森林图（按模块分面） ----------
cat("\n[2] 批量森林图\n")
keep_methods <- c("Inverse variance weighted", "Weighted median", "MR Egger")
plot_dat <- subset(res_all, method %in% keep_methods)
plot_dat$label <- paste0(plot_dat$trait, " | ", plot_dat$method)
p1 <- ggplot(plot_dat, aes(x = OR, y = reorder(label, OR), color = sig)) +
  geom_point(size = 3) + geom_errorbarh(aes(xmin = OR_lci, xmax = OR_uci), height = 0.2) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10() +
  scale_color_manual(values = c("P<0.05" = "darkred", "NS" = "grey50")) +
  labs(title = "MR 分析结果汇总森林图", x = "OR (log scale)", y = "") +
  theme_bw() + theme(axis.text.y = element_text(size = 8))
ggplot2::ggsave("04.figures/visual_forest_all.png", p1, width = 9, height = 6)
cat("    04.figures/visual_forest_all.png 已保存\n")

# ---------- 3. 多暴露热图（LDL/HDL/TG x 结局 OR 矩阵） ----------
cat("\n[3] 多暴露-多结局热图\n")
multi <- read.csv("02.analysis/opengwas/multi/comparison_summary.csv", stringsAsFactors = FALSE)
ivw_multi <- subset(multi, method == "Inverse variance weighted")
heat <- data.frame(
  暴露 = ivw_multi$trait,
  结局 = "CHD",
  OR = ivw_multi$OR,
  P = ivw_multi$pval)
heat$label <- sprintf("%.2f (P=%.1g)", heat$OR, heat$P)
p2 <- ggplot(heat, aes(x = 结局, y = 暴露, fill = log(OR))) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "log(OR)") +
  labs(title = "多暴露 MR 热图: LDL/HDL/TG -> CHD (IVW)", x = "", y = "") +
  theme_minimal()
ggplot2::ggsave("04.figures/visual_heatmap_multi.png", p2, width = 6, height = 4)
cat("    04.figures/visual_heatmap_multi.png 已保存\n")

# ---------- 4. MR 扫描热图（LDL x 多结局） ----------
cat("\n[4] MR 扫描热图（LDL-C x 多结局）\n")
scan <- read.csv("02.analysis/mr_scan/scan_results.csv", stringsAsFactors = FALSE)
ivw_scan <- subset(scan, method == "Inverse variance weighted")
heat2 <- data.frame(暴露 = "LDL-C", 结局 = ivw_scan$trait,
                    OR = ivw_scan$OR, P = ivw_scan$pval)
heat2$label <- sprintf("%.2f\n(P=%.1g)", heat2$OR, heat2$P)
p3 <- ggplot(heat2, aes(x = 结局, y = 暴露, fill = log(OR))) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3.5) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "log(OR)") +
  labs(title = "MR 扫描热图: LDL-C 对多结局因果效应 (IVW)", x = "", y = "") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggplot2::ggsave("04.figures/visual_heatmap_scan.png", p3, width = 8, height = 4)
cat("    04.figures/visual_heatmap_scan.png 已保存\n")

# ---------- 5. 汇总表 ----------
cat("\n[5] 可视化汇总表\n")
summary_tbl <- data.frame(
  图 = c("visual_forest_all.png", "visual_heatmap_multi.png", "visual_heatmap_scan.png"),
  内容 = c("各模块 MR 结果森林图（IVW/加权中位数/Egger）",
           "LDL/HDL/TG -> CHD 热图", "LDL-C -> 6 结局扫描热图"))
write.csv(summary_tbl, "02.analysis/visual/report_summary.csv", row.names = FALSE)
print(summary_tbl, row.names = FALSE)
cat("\n可视化报告完成 ✔\n")
