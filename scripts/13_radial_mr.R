#!/usr/bin/env Rscript
# 13_radial_mr.R — 径向 MR（Radial MR）：基于 Q 统计量的离群值检验
#
# 数据: 02.analysis/harmonised.csv（Telomere_length -> CHD，17 个工具变量）
# 方法: Radial IVW / Radial Egger（Bowden et al. 2018）
#   - 径向尺度上每个 SNP 贡献 Q 统计量，Q 超出卡方阈值者为离群 SNP
#   - 剔除离群 SNP 后重新估计，结果更稳健
#
# 输出: 02.analysis/radial/ 下结果与图

options(width = 150)
dir.create("02.analysis/radial", showWarnings = FALSE)
suppressMessages(library(RadialMR))

cat("========================================\n")
cat("径向 MR（Radial MR, Bowden et al. 2018）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

dat <- read.csv("02.analysis/harmonised.csv", stringsAsFactors = FALSE)
dat <- subset(dat, mr_keep)
cat("[1] 工具变量数:", nrow(dat), "\n")

# 转换为 RadialMR 格式（data.frame rmr_format）
r_input <- tsmr_to_rmr_format(dat)

# ---------- Radial IVW ----------
cat("\n[2] Radial IVW（含 Q 统计量离群检验）\n")
ivw <- ivw_radial(r_input, alpha = 0.05)
cat(sprintf("    Q 统计量总和 = %.2f,  df = %d,  P = %.3g\n",
            ivw$qstatistic, ivw$df, pchisq(ivw$qstatistic, ivw$df, lower.tail = FALSE)))
cat(sprintf("    估计效应(Mod.2nd) = %.4f (SE=%.4f, P=%.3g)\n",
            ivw$coef[1, 1], ivw$coef[1, 2], ivw$coef[1, 4]))
out_ivw <- ivw$outliers
cat("    离群 SNP:", ifelse(nrow(out_ivw) == 0, "无", paste(out_ivw$SNP, collapse = ", ")), "\n")

# 剔除离群后重估
if (nrow(out_ivw) > 0) {
  r_clean <- r_input[!(r_input$SNP %in% out_ivw$SNP), ]
  ivw2 <- ivw_radial(r_clean, alpha = 0.05)
  cat(sprintf("    剔除 %d 个离群 SNP 后: 效应=%.4f (SE=%.4f, P=%.3g)\n",
              nrow(out_ivw), ivw2$coef[1, 1], ivw2$coef[1, 2], ivw2$coef[1, 4]))
} else {
  cat("    无需剔除离群 SNP\n")
}

# ---------- Radial Egger ----------
cat("\n[3] Radial Egger（检验定向多效性）\n")
egg <- egger_radial(r_input, alpha = 0.05)
cat(sprintf("    Q 统计量 = %.2f,  df = %d,  P = %.3g\n",
            egg$qstatistic, egg$df, pchisq(egg$qstatistic, egg$df, lower.tail = FALSE)))
cat(sprintf("    截距(多效性) = %.4f (SE=%.4f, P=%.3g)\n",
            egg$coef[1, 1], egg$coef[1, 2], egg$coef[1, 4]))
cat(sprintf("    效应(斜率) = %.4f (SE=%.4f, P=%.3g)\n",
            egg$coef[2, 1], egg$coef[2, 2], egg$coef[2, 4]))

# ---------- 结果保存 ----------
res <- data.frame(
  方法 = c("Radial IVW(Mod.2nd)", "Radial Egger-截距(多效性)", "Radial Egger-效应"),
  估计 = c(ivw$coef[1, 1], egg$coef[1, 1], egg$coef[2, 1]),
  SE = c(ivw$coef[1, 2], egg$coef[1, 2], egg$coef[2, 2]),
  P = c(ivw$coef[1, 4], egg$coef[1, 4], egg$coef[2, 4]),
  Q = c(ivw$qstatistic, egg$qstatistic, egg$qstatistic),
  Q_P = c(pchisq(ivw$qstatistic, ivw$df, lower.tail = FALSE),
          pchisq(egg$qstatistic, egg$df, lower.tail = FALSE),
          pchisq(egg$qstatistic, egg$df, lower.tail = FALSE)))
write.csv(res, "02.analysis/radial/results.csv", row.names = FALSE)
if (nrow(out_ivw) > 0) write.csv(out_ivw, "02.analysis/radial/outliers.csv", row.names = FALSE)

# ---------- 径向图 ----------
cat("\n[4] 径向图\n")
p <- tryCatch(plot_radial(ivw, radial_scale = TRUE), error = function(e) {
  cat("    plotly 径向图失败:", conditionMessage(e), "\n"); NULL
})
if (!is.null(p)) {
  ok <- tryCatch({
    htmlwidgets::saveWidget(p, "02.analysis/radial/radial_ivw.html", selfcontained = TRUE)
    TRUE
  }, error = function(e) {
    cat("    HTML 保存失败（不影响分析结果）:", conditionMessage(e), "\n"); FALSE
  })
  if (ok) cat("    02.analysis/radial/radial_ivw.html 已保存\n")
}
cat("\n径向 MR 完成 ✔\n")
