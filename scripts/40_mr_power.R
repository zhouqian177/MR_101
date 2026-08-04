#!/usr/bin/env Rscript
# 40_mr_power.R — MR 统计功效计算（Power Calculation）
#
# 背景: MR 研究的功效（power）取决于: 工具变量强度(R^2/F)、样本量(N)、
#       真实效应大小(beta)、工具变量个数(K)。功效计算帮助研究者:
#       1) 评估现有数据能否检测到预期效应
#       2) 规划所需样本量/工具变量数
#       3) 解释阴性结果是否因功效不足
#
# 公式（Brion et al. 2013, mRnd）:
#   非中心参数 λ = N * K * F * beta^2 / (F * K + N)
#   或 λ = N * R^2 * beta^2 / (1 - R^2)  (单工具)
#   功效 = 1 - pchisq(χ²_crit, 1, λ)
#
# 输出: 02.analysis/power/ 表 + 04.figures/power_*.png

options(width = 150)
dir.create("02.analysis/power", showWarnings = FALSE)
set.seed(20260804)

cat("========================================\n")
cat("MR 统计功效计算\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 功效函数 ----------
mr_power <- function(N, K, R2, beta, alpha = 0.05) {
  # 非中心参数（Burgess 2014 公式）
  lambda <- N * K * R2 * beta^2 / (K * R2 * beta^2 + N)
  crit <- qchisq(alpha, 1, lower.tail = FALSE)
  1 - pchisq(crit, 1, lambda)
}

# ---------- 绘制功效曲线 ----------
cat("[1] 绘制功效曲线\n")

# 场景1: 固定样本量, 变化工具强度与效应大小
png("04.figures/power_curve_N.png", width = 800, height = 600, res = 110)
Ns <- c(5000, 10000, 50000, 100000, 200000)
K <- 30; R2 <- 0.03; betas <- seq(0.01, 0.30, 0.01)
plot(0, 0, type = "n", xlim = range(betas), ylim = c(0, 1),
     xlab = "真实因果效应 beta", ylab = "统计功效",
     main = paste0("MR 功效曲线 (K=", K, ", R²=", R2, ")"))
abline(h = 0.8, lty = 2, col = "grey")
for (n in Ns) {
  pw <- sapply(betas, function(b) mr_power(n, K, R2, b))
  lines(betas, pw, col = rainbow(length(Ns))[which(Ns == n)], lwd = 2)
  text(max(betas), tail(pw, 1), paste0("N=", n), cex = 0.7, pos = 4)
}
legend("topleft", paste0("N=", Ns), col = rainbow(length(Ns)), lwd = 2, cex = 0.7)
dev.off()
cat("    04.figures/power_curve_N.png 已保存\n")

# 场景2: 固定效应, 变化工具变量数(F)
png("04.figures/power_curve_K.png", width = 800, height = 600, res = 110)
beta_true <- 0.1; N <- 100000
Ks <- seq(5, 100, 5)
R2s <- c(0.01, 0.02, 0.05, 0.10)
plot(0, 0, type = "n", xlim = range(Ks), ylim = c(0, 1),
     xlab = "工具变量数 K", ylab = "统计功效",
     main = paste0("MR 功效曲线 (beta=", beta_true, ", N=", N, ")"))
abline(h = 0.8, lty = 2, col = "grey")
for (r in R2s) {
  pw <- sapply(Ks, function(k) mr_power(N, k, r, beta_true))
  lines(Ks, pw, col = rainbow(length(R2s))[which(R2s == r)], lwd = 2)
  text(max(Ks), tail(pw, 1), paste0("R²=", r), cex = 0.7, pos = 4)
}
legend("topleft", paste0("R²=", R2s), col = rainbow(length(R2s)), lwd = 2, cex = 0.7)
dev.off()
cat("    04.figures/power_curve_K.png 已保存\n")

# ---------- 实际数据示例 ----------
cat("\n[2] 实际数据示例: LDL-C -> CHD 功效\n")
# 从在线分析结果中获取参数
d <- read.csv("02.analysis/opengwas/online/mr_results.csv", stringsAsFactors = FALSE)
ivw <- subset(d, method == "Inverse variance weighted")
N_eff <- 60801 * 123504 / (60801 + 123504)  # CHD 有效样本量
K_use <- ivw$nsnp; beta_use <- abs(ivw$b)
# 估计 R²（从 F 统计量近似: R² = F*K/(N + F*K)）
F_mean <- 177  # 之前 LDL 数据 F 均值
R2_est <- F_mean * K_use / (N_eff + F_mean * K_use)
power_actual <- mr_power(N_eff, K_use, R2_est, beta_use)
cat(sprintf("    LDL-C -> CHD 实际功效: %.1f%% (K=%d, R²=%.3f, beta=%.3f, N=%.0f)\n",
            power_actual * 100, K_use, R2_est, beta_use, N_eff))

# 所需最小样本量（达到 80% 功效）
cat(sprintf("    需要 %.0f 样本达到 80%% 功效（当前 %.0f）\n",
            N_eff / power_actual * 0.8, N_eff))

# ---------- 保存结果 ----------
res <- data.frame(
  参数 = c("K (工具变量数)", "N (有效样本量)", "R² (暴露解释方差)",
           "beta (真实效应)", "统计功效"),
  值 = c(K_use, round(N_eff), round(R2_est, 4), round(beta_use, 4), round(power_actual, 3)))
write.csv(res, "02.analysis/power/power_ldl_chd.csv", row.names = FALSE)
print(res, row.names = FALSE)
cat("\nMR 功效计算完成 ✔\n")