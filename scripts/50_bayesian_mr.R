#!/usr/bin/env Rscript
# 50_bayesian_mr.R — 贝叶斯模型平均 MR（Bayesian model averaging）
#
# 实现一个简单的贝叶斯 MR：对每个工具变量的 Wald ratio 做贝叶斯元分析，
# 用网格近似计算后验分布，自动权衡有效/无效 IV 的贡献。
#
# 方法:
#   1) 每个工具变量计算 Wald ratio theta_i = beta_Y / beta_X
#   2) 贝叶斯模型: theta_i ~ N(theta, se_i^2 + tau^2) 随机效应
#   3) 先验: theta ~ N(0, 10), tau ~ Half-Cauchy(0, 1)
#   4) 用网格近似计算后验均值与 95% 可信区间
#
# 数据: LDL-C -> CHD（28 个工具变量，在线数据）
# 输出: 02.analysis/bayesian_mr/ 结果

options(width = 150)
dir.create("02.analysis/bayesian_mr", showWarnings = FALSE)
set.seed(20260804)

cat("========================================\n")
cat("贝叶斯模型平均 MR\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 数据 ----------
cat("[1] 数据: LDL-C -> CHD（28 个工具变量）\n")
d <- read.csv("02.analysis/opengwas/online/harmonised.csv", stringsAsFactors = FALSE)
d <- subset(d, mr_keep)
theta_i <- d$beta.outcome / d$beta.exposure
se_i <- d$se.outcome / abs(d$beta.exposure)
cat("    工具变量数:", length(theta_i), "\n")
cat("    Wald ratio 范围:", sprintf("%.3f ~ %.3f", min(theta_i), max(theta_i)), "\n")

# ---------- 贝叶斯元分析（网格近似）----------
cat("\n[2] 贝叶斯元分析（网格近似）\n")
theta_grid <- seq(-2, 2, length.out = 2000)
prior_theta <- dnorm(theta_grid, 0, 10)  # 弱先验 N(0, 10)
tau_grid <- seq(0.01, 1, length.out = 500)
prior_tau <- 2 / (pi * (1 + tau_grid^2))  # Half-Cauchy(0,1)

# 计算联合后验（网格近似）
log_lik <- function(theta, tau) {
  sum(dnorm(theta_i, theta, sqrt(se_i^2 + tau^2), log = TRUE))
}
post <- matrix(0, length(theta_grid), length(tau_grid))
for (i in seq_along(theta_grid)) {
  for (j in seq_along(tau_grid)) {
    post[i, j] <- log_lik(theta_grid[i], tau_grid[j]) +
      log(prior_theta[i]) + log(prior_tau[j])
  }
}
post <- exp(post - max(post))
post <- post / sum(post)

# 边缘后验
post_theta <- rowSums(post)
post_theta <- post_theta / sum(post_theta)

# 后验均值与可信区间
cdf <- cumsum(post_theta)
bayes_mean <- sum(theta_grid * post_theta)
bayes_lower <- theta_grid[which.min(abs(cdf - 0.025))]
bayes_upper <- theta_grid[which.min(abs(cdf - 0.975))]

cat(sprintf("    后验均值: %.4f\n", bayes_mean))
cat(sprintf("    95%% 可信区间: (%.4f, %.4f)\n", bayes_lower, bayes_upper))
cat(sprintf("    后验 P(theta>0): %.3f\n", sum(post_theta[theta_grid > 0])))

# ---------- 与 IVW 对比 ----------
cat("\n[3] 与经典方法对比\n")
w <- 1 / se_i^2
ivw_b <- sum(theta_i * w) / sum(w)
ivw_se <- sqrt(1 / sum(w))
cat(sprintf("    IVW: beta=%.4f (SE=%.4f)\n", ivw_b, ivw_se))
cat(sprintf("    贝叶斯: beta=%.4f (95%%CI: %.4f-%.4f)\n", bayes_mean, bayes_lower, bayes_upper))

# ---------- 后验分布图 ----------
cat("\n[4] 绘图: 后验分布\n")
png("04.figures/bayesian_posterior.png", width = 800, height = 600, res = 110)
plot(theta_grid, post_theta, type = "l", lwd = 2, col = "darkred",
     xlab = "因果效应 theta", ylab = "后验概率密度",
     main = "贝叶斯 MR: 后验分布 (LDL-C -> CHD)")
abline(v = bayes_mean, col = "darkred", lwd = 2, lty = 2)
abline(v = bayes_lower, col = "grey", lty = 3)
abline(v = bayes_upper, col = "grey", lty = 3)
legend("topright", c("后验均值", "95% 可信区间"),
       col = c("darkred", "grey"), lty = c(2, 3), lwd = 2, bty = "n")
dev.off()
cat("    04.figures/bayesian_posterior.png 已保存\n")

res <- data.frame(
  方法 = c("IVW", "Bayesian MR"),
  beta = c(ivw_b, bayes_mean),
  lower = c(ivw_b - 1.96 * ivw_se, bayes_lower),
  upper = c(ivw_b + 1.96 * ivw_se, bayes_upper))
write.csv(res, "02.analysis/bayesian_mr/results.csv", row.names = FALSE)
print(res, row.names = FALSE)
cat("\n贝叶斯模型平均 MR 完成 ✔\n")