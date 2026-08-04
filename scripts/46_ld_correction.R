#!/usr/bin/env Rscript
# 46_ld_correction.R — MR 的 LD 矩阵校正
#
# 背景: 标准两样本 MR 假设工具变量间相互独立（无 LD）。
#       当工具变量存在 LD 时，IVW 估计会因协方差未校正而产生偏倚。
#       MendelianRandomization::mr_ivw 的 correl=TRUE 参数可接受
#       LD 矩阵（SNP 间相关系数矩阵）进行校正。
#
# 演示: 模拟有 LD 的 SNP 数据，比较未校正与校正后的 IVW 估计
#
# 输出: 02.analysis/ld_correction/ 结果

options(width = 150)
dir.create("02.analysis/ld_correction", showWarnings = FALSE)
set.seed(20260804)
suppressMessages(library(MendelianRandomization))

cat("========================================\n")
cat("MR 的 LD 矩阵校正\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 模拟有 LD 的数据 ----------
cat("[1] 模拟数据（30 个工具变量，其中 10 对有 LD，r=0.6）\n")
n <- 50000; K <- 30
# 先模拟独立多变量正态分布的基因型（部分 SNP 有 LD）
library(MASS)
Sigma <- diag(K)
# 设置 10 对 SNP 有 LD（r=0.6）
for (i in 1:10) {
  Sigma[i*2-1, i*2] <- 0.6; Sigma[i*2, i*2-1] <- 0.6
}
G <- mvrnorm(n, mu = rep(0, K), Sigma = Sigma)
G <- pmin(pmax(round(G + 2), 0), 2)  # 离散化为 0/1/2 剂量
beta_g <- runif(K, 0.05, 0.10)
X <- as.vector(G %*% beta_g) + rnorm(n)
Y <- 0.3 * X + rnorm(n)

# 汇总统计量
beta_x <- se_x <- beta_y <- se_y <- numeric(K)
for (k in 1:K) {
  f1 <- summary(lm(X ~ G[, k])); f2 <- summary(lm(Y ~ G[, k]))
  beta_x[k] <- f1$coef[2, 1]; se_x[k] <- f1$coef[2, 2]
  beta_y[k] <- f2$coef[2, 1]; se_y[k] <- f2$coef[2, 2]
}

# LD 矩阵（SNP 间相关系数）
ld_matrix <- cor(G)
cat("    LD 矩阵维度:", nrow(ld_matrix), "x", ncol(ld_matrix), "\n")
cat("    LD 矩阵中 |r|>0.3 的比例:", 
    round(sum(abs(ld_matrix[lower.tri(ld_matrix)]) > 0.3) / sum(lower.tri(ld_matrix)) * 100, 1), "%\n")

# ---------- 未校正 vs 校正 ----------
cat("\n[2] IVW 未校正 LD（standard）\n")
mrin <- mr_input(bx = beta_x, bxse = se_x, by = beta_y, byse = se_y)
ivw_uncor <- mr_ivw(mrin)
cat(sprintf("    beta=%.4f (SE=%.4f, P=%.3g)  [真实 0.3]\n",
            ivw_uncor@Estimate, ivw_uncor@StdError, ivw_uncor@Pvalue))

cat("\n[3] IVW 校正 LD（correl=TRUE, mr_input 传入 correlation）\n")
mrin_ld <- mr_input(bx = beta_x, bxse = se_x, by = beta_y, byse = se_y,
                    correlation = ld_matrix)
ivw_cor <- mr_ivw(mrin_ld, correl = TRUE)
cat(sprintf("    beta=%.4f (SE=%.4f, P=%.3g)  [真实 0.3]\n",
            ivw_cor@Estimate, ivw_cor@StdError, ivw_cor@Pvalue))

# ---------- 对比 ----------
cat("\n[4] 对比\n")
cat(sprintf("    未校正: SE=%.4f\n", ivw_uncor@StdError))
cat(sprintf("    校正后: SE=%.4f\n", ivw_cor@StdError))
cat(sprintf("    SE 变化: %.1f%%\n",
            (ivw_cor@StdError / ivw_uncor@StdError - 1) * 100))
cat("    结论: LD 校正影响标准误（SE），校正后更准确\n")

res <- data.frame(
  方法 = c("IVW (未校正 LD)", "IVW (校正 LD)"),
  beta = c(ivw_uncor@Estimate, ivw_cor@Estimate),
  SE = c(ivw_uncor@StdError, ivw_cor@StdError),
  P = c(ivw_uncor@Pvalue, ivw_cor@Pvalue))
write.csv(res, "02.analysis/ld_correction/results.csv", row.names = FALSE)
print(res, row.names = FALSE)
cat("\nLD 矩阵校正完成 ✔\n")