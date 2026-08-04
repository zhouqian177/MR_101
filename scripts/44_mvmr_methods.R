#!/usr/bin/env Rscript
# 44_mvmr_methods.R — MR 多变量方法补充（MVMedian/MVMaxLik/MRMaxLik）
#
# 补充 MendelianRandomization 包中尚未演示的方法:
#   1) mr_mvmedian — 多变量加权中位数（对离群稳健）
#   2) mr_maxlik  — 最大似然 MR（基于似然比，效率高）
#   3) mr_mvivw   — 多变量 IVW（已有 11_mvmr.R, 此处对比参考）
#
# 数据: 模拟双暴露（X1, X2）→ 结局 Y
# 真实效应: X1→Y=0.4, X2→Y=-0.3
# 输出: 02.analysis/mvmr_supp/ 结果

options(width = 150)
dir.create("02.analysis/mvmr_supp", showWarnings = FALSE)
set.seed(20260804)
suppressMessages(library(MendelianRandomization))

cat("========================================\n")
cat("MR 多变量方法补充\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 模拟数据 ----------
cat("[1] 模拟数据（X1→Y=0.4, X2→Y=-0.3, 30 工具变量）\n")
n <- 30000; K1 <- 15; K2 <- 15
G1 <- matrix(rbinom(n * K1, 2, 0.3), n, K1)
G2 <- matrix(rbinom(n * K2, 2, 0.3), n, K2)
b1 <- runif(K1, 0.05, 0.10); b2 <- runif(K2, 0.05, 0.10)
U <- rnorm(n)
X1 <- as.vector(G1 %*% b1) + 0.5 * U + rnorm(n)
X2 <- as.vector(G2 %*% b2) + 0.5 * U + rnorm(n)
Y  <- 0.4 * X1 - 0.3 * X2 + rnorm(n)

# 汇总统计量
K <- K1 + K2
beta_x1 <- beta_x2 <- beta_y <- se_x1 <- se_x2 <- se_y <- numeric(K)
for (k in 1:K1) {
  g <- G1[, k]
  beta_x1[k] <- coef(summary(lm(X1 ~ g)))[2,1]; se_x1[k] <- coef(summary(lm(X1 ~ g)))[2,2]
  beta_x2[k] <- coef(summary(lm(X2 ~ g)))[2,1]; se_x2[k] <- coef(summary(lm(X2 ~ g)))[2,2]
  beta_y[k]  <- coef(summary(lm(Y  ~ g)))[2,1]; se_y[k]  <- coef(summary(lm(Y  ~ g)))[2,2]
}
for (k in 1:K2) {
  g <- G2[, k]
  beta_x1[K1+k] <- coef(summary(lm(X1 ~ g)))[2,1]; se_x1[K1+k] <- coef(summary(lm(X1 ~ g)))[2,2]
  beta_x2[K1+k] <- coef(summary(lm(X2 ~ g)))[2,1]; se_x2[K1+k] <- coef(summary(lm(X2 ~ g)))[2,2]
  beta_y[K1+k]  <- coef(summary(lm(Y  ~ g)))[2,1]; se_y[K1+k]  <- coef(summary(lm(Y  ~ g)))[2,2]
}

mvin <- mr_mvinput(bx = cbind(beta_x1, beta_x2),
                   bxse = cbind(se_x1, se_x2),
                   by = beta_y, byse = se_y,
                   exposure = c("X1", "X2"), outcome = "Y")

# ---------- 各方法 ----------
cat("\n[2] MV-IVW（多变量逆方差加权，参考）\n")
ivw <- mr_mvivw(mvin)
cat(sprintf("    X1: %.4f (SE=%.4f, P=%.3g)  X2: %.4f (SE=%.4f, P=%.3g)\n",
            ivw@Estimate[1], ivw@StdError[1], ivw@Pvalue[1],
            ivw@Estimate[2], ivw@StdError[2], ivw@Pvalue[2]))

cat("\n[3] MV-Median（多变量加权中位数，对离群稳健）\n")
med <- mr_mvmedian(mvin)
cat(sprintf("    X1: %.4f (SE=%.4f, P=%.3g)  X2: %.4f (SE=%.4f, P=%.3g)\n",
            med@Estimate[1], med@StdError[1], med@Pvalue[1],
            med@Estimate[2], med@StdError[2], med@Pvalue[2]))

cat("\n[4] 单变量 MR 对比（MaxLik - 最大似然法）\n")
mrin_x1 <- mr_input(bx = beta_x1, bxse = se_x1, by = beta_y, byse = se_y)
ml1 <- mr_maxlik(mrin_x1)
cat(sprintf("    X1 MaxLik: beta=%.4f (SE=%.4f, P=%.3g)  [真实 0.4]\n",
            ml1@Estimate, ml1@StdError, ml1@Pvalue))

# ---------- 汇总 ----------
cat("\n[5] 汇总对比\n")
res <- data.frame(
  方法 = c("MV-IVW", "MV-IVW", "MV-Median", "MV-Median", "MaxLik-X1"),
  暴露 = c("X1", "X2", "X1", "X2", "X1"),
  估计 = c(ivw@Estimate[1], ivw@Estimate[2], med@Estimate[1], med@Estimate[2], ml1@Estimate),
  真实值 = c(0.4, -0.3, 0.4, -0.3, 0.4))
print(res, row.names = FALSE)
write.csv(res, "02.analysis/mvmr_supp/results.csv", row.names = FALSE)
cat("\nMR 多变量方法补充完成 ✔\n")