#!/usr/bin/env Rscript
# 15_coloc.R — 共定位分析（Colocalization）：区分"共享因果变异"与"连锁不平衡"
#
# 背景: MR 发现的暴露-结局关联可能源于同一因果变异（共定位）或两个紧密连锁
# 的不同变异（LD 假关联）。coloc 通过贝叶斯方法计算五种假设的后验概率:
#   H0: 无关联  H1: 仅暴露关联  H2: 仅结局关联  H3: 两关联但变异不同  H4: 共享因果变异
# PP.H4 高（通常 > 0.8）支持共定位。
#
# 数据: coloc 包内置 coloc_test_data（D1~D4 为不同模拟区域）
#
# 输出: 02.analysis/coloc/ 下结果

options(width = 150)
dir.create("02.analysis/coloc", showWarnings = FALSE)
suppressMessages(library(coloc))

cat("========================================\n")
cat("共定位分析（coloc.abf）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

data("coloc_test_data")
cat("[1] 内置数据区域: D1~D4，每区域含 beta/varbeta/snp/MAF/N 等\n")

run_coloc <- function(d1, d2, label) {
  cat("\n----------------------------------------\n")
  cat("[2] 区域:", label, "\n")
  res <- coloc.abf(d1, d2)
  pp <- res$summary["PP.H4.abf"]
  # 输出关键后验概率
  cat(sprintf("    PP.H0=%.3f PP.H1=%.3f PP.H2=%.3f PP.H3=%.3f PP.H4=%.3f\n",
              res$summary["PP.H0.abf"], res$summary["PP.H1.abf"],
              res$summary["PP.H2.abf"], res$summary["PP.H3.abf"],
              res$summary["PP.H4.abf"]))
  concl <- ifelse(pp > 0.8, "支持共定位（共享因果变异）",
           ifelse(pp > 0.5, "中等证据", "不支持共定位（可能为 LD 假关联）"))
  cat("    结论:", concl, "\n")
  list(res = res, label = label, concl = concl)
}

# D1 vs D2: 已知模拟为共定位区域
r1 <- run_coloc(coloc_test_data$D1, coloc_test_data$D2, "D1 vs D2（模拟: 共定位）")
# D1 vs D3 / D1 vs D4: 对比区域
r2 <- run_coloc(coloc_test_data$D1, coloc_test_data$D3, "D1 vs D3（模拟: 无共定位）")
r3 <- run_coloc(coloc_test_data$D1, coloc_test_data$D4, "D1 vs D4（模拟: 部分信号）")

# ---------- 结果汇总 ----------
res <- data.frame(
  区域 = c(r1$label, r2$label, r3$label),
  PP.H3 = c(r1$res$summary["PP.H3.abf"], r2$res$summary["PP.H3.abf"], r3$res$summary["PP.H3.abf"]),
  PP.H4 = c(r1$res$summary["PP.H4.abf"], r2$res$summary["PP.H4.abf"], r3$res$summary["PP.H4.abf"]),
  结论 = c(r1$concl, r2$concl, r3$concl))
print(res, row.names = FALSE)
write.csv(res, "02.analysis/coloc/results.csv", row.names = FALSE)

# ---------- 敏感性分析（先验 p12 变化） ----------
cat("\n[3] 敏感性分析: p12 先验从 1e-6 到 1e-4\n")
sens <- sensitivity(r1$res, rule = "H4 > 0.8")
cat("    PP.H4 随 p12 变化范围:",
    sprintf("%.3f ~ %.3f", min(sens$PP.H4), max(sens$PP.H4)), "\n")
cat("\n共定位分析完成 ✔\n")
