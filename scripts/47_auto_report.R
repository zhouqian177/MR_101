#!/usr/bin/env Rscript
# 47_auto_report.R — MR 结果自动报告（论文格式汇总表 + 全景图）
#
# 自动汇总本仓库所有关键 MR 分析结果，生成论文格式的汇总表与全景图。
# 输出: 02.analysis/auto_report/ 汇总表 + 04.figures/report_*.png

options(width = 150)
dir.create("02.analysis/auto_report", showWarnings = FALSE)
suppressMessages(library(ggplot2))

cat("========================================\n")
cat("MR 结果自动报告\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 汇总各模块结果 ----------
cat("[1] 汇总各模块结果\n")
collect <- function(f, label) {
  if (!file.exists(f)) return(NULL)
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(d) == 0) return(NULL)
  if ("method" %in% colnames(d)) {
    sub <- subset(d, method %in% c("Inverse variance weighted", "Weighted median"))
    if (nrow(sub) == 0) return(NULL)
    sub <- sub[, intersect(colnames(sub), c("method", "b", "se", "pval", "OR", "OR_lci", "OR_uci", "nsnp"))]
    sub$分析 <- label
    sub$beta <- sub$b
    sub$P <- sub$pval
    sub$工具数 <- sub$nsnp
    return(sub[, c("分析", "method", "beta", "se", "P", "OR", "OR_lci", "OR_uci", "工具数")])
  }
  NULL
}
parts <- list(
  collect("02.analysis/mr_results.csv", "Telomere→CHD"),
  collect("02.analysis/opengwas/online/mr_results.csv", "LDL-C→CHD(在线)"),
  collect("02.analysis/mr_scan/scan_results.csv", "LDL-C→多结局(MR扫描)"),
  collect("02.analysis/opengwas/multi/comparison_summary.csv", "LDL/HDL/TG→CHD(多暴露)"),
  collect("02.analysis/panorama/panorama_results.csv", "方法全景对比"))
parts <- parts[!sapply(parts, is.null)]
res_all <- do.call(rbind, parts)
cat("    汇总行数:", nrow(res_all), "\n")

# ---------- 论文格式表格 ----------
cat("\n[2] 生成论文格式表格\n")
res_all$OR_str <- sprintf("%.2f (%.2f-%.2f)", res_all$OR, res_all$OR_lci, res_all$OR_uci)
res_all$P_str <- ifelse(res_all$P < 0.001, "<0.001",
                        sprintf("%.3f", res_all$P))
res_all$beta_str <- sprintf("%.3f", res_all$beta)
res_all$se_str <- sprintf("%.3f", res_all$se)
table_out <- res_all[, c("分析", "method", "工具数", "beta_str", "se_str", "OR_str", "P_str")]
colnames(table_out) <- c("分析", "方法", "n", "beta", "SE", "OR(95%CI)", "P")
write.csv(table_out, "02.analysis/auto_report/paper_ready_table.csv", row.names = FALSE)
cat("    已保存 02.analysis/auto_report/paper_ready_table.csv\n")
print(table_out, row.names = FALSE)

# ---------- 全景图 ----------
cat("\n[3] 全景图\n")
ivw_dat <- subset(res_all, method == "Inverse variance weighted" | method == "Weighted median")
ivw_dat$label <- paste0(ivw_dat$分析, " (", ivw_dat$method, ", n=", ivw_dat$工具数, ")")
p <- ggplot(ivw_dat, aes(x = OR, y = reorder(label, OR),
                          color = ifelse(P < 0.05, "P<0.05", "NS"))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = OR_lci, xmax = OR_uci), height = 0.2) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  scale_x_log10() +
  scale_color_manual(values = c("P<0.05" = "darkred", "NS" = "grey50")) +
  labs(title = "MR 仓库结果全景图: IVW 与加权中位数",
       subtitle = paste0("汇总 ", nrow(ivw_dat), " 个分析"),
       x = "OR (95% CI, log scale)", y = NULL, color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggplot2::ggsave("04.figures/report_overview.png", p, width = 10, height = 6)
cat("    04.figures/report_overview.png 已保存\n")

# ---------- 汇总统计 ----------
cat("\n[4] 汇总统计\n")
cat("    仓库分析脚本总数:", length(list.files("scripts", pattern = "\\.(R|py)$")), "\n")
cat("    仓库文档总数:", length(list.files("docs", pattern = "\\.md$")), "\n")
cat("    本次汇总分析数:", nrow(res_all), "\n")
cat("    显著结果数(IVW):", sum(ivw_dat$P < 0.05 & ivw_dat$method == "Inverse variance weighted"), "\n")
cat("\nMR 结果自动报告完成 ✔\n")