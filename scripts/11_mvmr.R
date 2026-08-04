#!/usr/bin/env Rscript
# 11_mvmr.R — 多变量 MR（MVMR）：同时估计多个暴露的因果效应
#
# 模拟数据:
#   X1 = G1*b1 + e1          （暴露1，15 个工具 SNP）
#   X2 = G2*b2 + e2          （暴露2，15 个工具 SNP，与 X1 相关 rho=0.3）
#   Y  = 0.4*X1 - 0.3*X2 + e （真实因果效应: X1→Y=0.4, X2→Y=-0.3）
#
# 演示:
#   1) 单变量 IVW（未校正 X2，估计 X1 时有偏）
#   2) 多变量 IVW / Egger / 中位数（校正暴露间相关性）
#
# 输出: 02.analysis/mvmr/ 下结果

options(width = 150)
dir.create("02.analysis/mvmr", showWarnings = FALSE)
set.seed(20260804)
suppressMessages(library(MendelianRandomization))

cat("========================================\n")
cat("多变量 MR（MVMR）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 模拟个体水平数据，生成汇总统计量 ----------
cat("[1] 模拟数据（X1→Y=0.4, X2→Y=-0.3, X1-X2 相关）\n")
n <- 30000                    # 个体数
K1 <- 15; K2 <- 15            # 每个暴露的工具 SNP 数
K <- K1 + K2
G1 <- matrix(rbinom(n * K1, 2, 0.3), n, K1)   # 暴露1 工具基因型
G2 <- matrix(rbinom(n * K2, 2, 0.3), n, K2)   # 暴露2 工具基因型
b1 <- runif(K1, 0.05, 0.10); b2 <- runif(K2, 0.05, 0.10)
U  <- rnorm(n)                                 # 共享混杂（使 X1、X2 相关）
X1 <- as.vector(G1 %*% b1) + 0.5 * U + rnorm(n)
X2 <- as.vector(G2 %*% b2) + 0.5 * U + rnorm(n)
Y  <- 0.4 * X1 - 0.3 * X2 + rnorm(n)
cat("    真实效应: X1→Y=0.4, X2→Y=-0.3; 工具 SNP:", K, "个\n")

# 每个 SNP 分别对 X1、X2、Y 回归，得到汇总统计量（beta/se）
beta_x1 <- beta_x2 <- beta_y <- se_x1 <- se_x2 <- se_y <- numeric(K)
for (k in 1:K1) {
  g <- G1[, k]
  f1 <- summary(lm(X1 ~ g)); f2 <- summary(lm(X2 ~ g)); fy <- summary(lm(Y ~ g))
  beta_x1[k] <- f1$coef[2, 1]; se_x1[k] <- f1$coef[2, 2]
  beta_x2[k] <- f2$coef[2, 1]; se_x2[k] <- f2$coef[2, 2]
  beta_y[k]  <- fy$coef[2, 1]; se_y[k]  <- fy$coef[2, 2]
}
for (k in 1:K2) {
  g <- G2[, k]
  f1 <- summary(lm(X1 ~ g)); f2 <- summary(lm(X2 ~ g)); fy <- summary(lm(Y ~ g))
  beta_x1[K1 + k] <- f1$coef[2, 1]; se_x1[K1 + k] <- f1$coef[2, 2]
  beta_x2[K1 + k] <- f2$coef[2, 1]; se_x2[K1 + k] <- f2$coef[2, 2]
  beta_y[K1 + k]  <- fy$coef[2, 1]; se_y[K1 + k]  <- fy$coef[2, 2]
}
cat("    汇总统计量生成完毕（", K, " 个 SNP × 3 个性状）\n\n")

# ---------- 2. 单变量 IVW（未校正，有偏演示） ----------
cat("[2] 单变量 IVW（忽略 X2，X1 估计应有偏）\n")
bw <- sum(beta_x1 * beta_y / se_y^2) / sum(beta_x1^2 / se_y^2)
cat(sprintf("    X1 单变量 IVW = %.3f （真实 0.4）<- 受 X2 混杂污染\n", bw))

# ---------- 3. 多变量 MR ----------
cat("\n[3] 多变量 MR（MendelianRandomization）\n")
mvin <- mr_mvinput(bx = cbind(beta_x1, beta_x2),
                   bxse = cbind(se_x1, se_x2),
                   by = beta_y, byse = se_y,
                   exposure = c("X1", "X2"), outcome = "Y",
                   snps = paste0("snp", 1:K))
res_ivw <- mr_mvivw(mvin)
res_egg <- mr_mvegger(mvin)
res_med <- mr_mvmedian(mvin)
cat("\n    [MV-IVW] 多变量逆方差加权\n")
print(res_ivw@Estimate, row.names = FALSE)
cat("\n    [MV-Egger] 多变量 Egger\n")
print(res_egg@Estimate, row.names = FALSE)
cat("\n    [MV-Median] 多变量加权中位数\n")
print(res_med@Estimate, row.names = FALSE)

# ---------- 4. 结果汇总对比 ----------
cat("\n[4] 结果对比（真实: X1=0.4, X2=-0.3）\n")
res <- data.frame(
  暴露 = c("X1", "X2"),
  真实效应 = c(0.4, -0.3),
  单变量IVW = c(bw, NA),
  MV_IVW = c(res_ivw@Estimate[1], res_ivw@Estimate[2]),
  MV_IVW_SE = c(res_ivw@StdError[1], res_ivw@StdError[2]))
print(res, row.names = FALSE)
write.csv(res, "02.analysis/mvmr/results.csv", row.names = FALSE)
cat("\n多变量 MR 完成 ✔\n")
