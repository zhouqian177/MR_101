#!/usr/bin/env Rscript
# 33_mvmr_mediation.R — 多变量中介 MR（MVMR-based mediation）
#
# 思路（Carter et al. 2021 等框架）:
#   暴露 X -> 中介 M -> 结局 Y, 且 X 对 Y 有直接效应
#   总效应 = 直接效应(X→Y, 校正M) + 间接效应(X→M→Y)
#   - 总效应:   用 X 工具变量做单变量 MR(X→Y)
#   - 直接效应: MVMR 同时纳入 X 与 M 的工具变量, 得 X 的校正系数
#   - 间接效应 = 总效应 - 直接效应（差值法）
#   - 间接效应也可 = a(X→M 用X工具) × b(M→Y 校正X后)
#
# 模拟设计（真实参数）:
#   X = Gx*bx + 0.5U + e1          暴露（15 工具）
#   M = 0.6X + Gm*bm + 0.5U + e2   中介（a=0.6, M 有 15 个工具）
#   Y = 0.3X + 0.5M + 0.5U + e3    结局（直接 c'=0.3, b=0.5）
#   真实: 总效应=0.3+0.6*0.5=0.60, 直接=0.30, 间接=0.30, 中介比例=50%
#
# 输出: 02.analysis/mvmr_mediation/ 结果

options(width = 150)
dir.create("02.analysis/mvmr_mediation", showWarnings = FALSE)
set.seed(20260804)
suppressMessages(library(MendelianRandomization))

cat("========================================\n")
cat("多变量中介 MR（MVMR-based Mediation）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 模拟个体数据 ----------
cat("[1] 模拟数据（a=0.6, b=0.5, c'=0.3, 中介比例=50%）\n")
n <- 30000; Kx <- 15; Km <- 15
Gx <- matrix(rbinom(n * Kx, 2, 0.3), n, Kx)
Gm <- matrix(rbinom(n * Km, 2, 0.3), n, Km)
bx <- runif(Kx, 0.06, 0.12); bm <- runif(Km, 0.06, 0.12)
U <- rnorm(n)
X <- as.vector(Gx %*% bx) + 0.5 * U + rnorm(n)
M <- 0.6 * X + as.vector(Gm %*% bm) + 0.5 * U + rnorm(n)
Y <- 0.3 * X + 0.5 * M + 0.5 * U + rnorm(n)
cat("    个体数:", n, "| 工具: X 用", Kx, "SNP, M 用", Km, "SNP\n")

# ---------- 2. 生成汇总统计量 ----------
cat("[2] 生成汇总统计量\n")
K <- Kx + Km
beta_x <- beta_m <- beta_y <- se_x <- se_m <- se_y <- numeric(K)
for (k in 1:Kx) {
  g <- Gx[, k]
  beta_x[k] <- coef(summary(lm(X ~ g)))[2, 1]; se_x[k] <- coef(summary(lm(X ~ g)))[2, 2]
  beta_m[k] <- coef(summary(lm(M ~ g)))[2, 1]; se_m[k] <- coef(summary(lm(M ~ g)))[2, 2]
  beta_y[k] <- coef(summary(lm(Y ~ g)))[2, 1]; se_y[k] <- coef(summary(lm(Y ~ g)))[2, 2]
}
for (k in 1:Km) {
  g <- Gm[, k]; j <- Kx + k
  beta_x[j] <- coef(summary(lm(X ~ g)))[2, 1]; se_x[j] <- coef(summary(lm(X ~ g)))[2, 2]
  beta_m[j] <- coef(summary(lm(M ~ g)))[2, 1]; se_m[j] <- coef(summary(lm(M ~ g)))[2, 2]
  beta_y[j] <- coef(summary(lm(Y ~ g)))[2, 1]; se_y[j] <- coef(summary(lm(Y ~ g)))[2, 2]
}
cat("    汇总统计量:", K, "SNP × 3 性状\n")

# ---------- 3. 总效应（单变量 MR: X→Y, 用 X 工具） ----------
cat("\n[3] 总效应: MR(X -> Y) 用 X 工具\n")
ivw <- function(bx, by, sey) {
  w <- 1 / sey^2
  b <- sum(bx * by * w) / sum(bx^2 * w); se <- sqrt(1 / sum(bx^2 * w))
  c(est = b, se = se, p = 2 * pnorm(-abs(b / se)))
}
total <- ivw(beta_x[1:Kx], beta_y[1:Kx], se_y[1:Kx])
cat(sprintf("    总效应 = %.3f (SE=%.3f)  [真实 0.60]\n", total["est"], total["se"]))

# ---------- 4. 直接效应（MVMR: X+M → Y, 校正后 X 系数） ----------
cat("\n[4] 直接效应: MVMR(X+M -> Y)\n")
mvin <- mr_mvinput(bx = cbind(beta_x, beta_m),
                   bxse = cbind(se_x, se_m),
                   by = beta_y, byse = se_y,
                   exposure = c("X", "M"), outcome = "Y")
res_mv <- mr_mvivw(mvin)
direct <- res_mv@Estimate[1]           # X 的校正系数 = 直接效应
cat(sprintf("    MVMR 校正后 X 效应(直接) = %.3f (SE=%.3f)  [真实 0.30]\n",
            direct, res_mv@StdError[1]))
cat(sprintf("    MVMR 校正后 M 效应(b')    = %.3f (SE=%.3f)  [真实 0.50]\n",
            res_mv@Estimate[2], res_mv@StdError[2]))

# ---------- 5. 间接效应与中介比例 ----------
cat("\n[5] 中介效应分解（差值法）\n")
indirect <- total["est"] - direct
prop <- indirect / total["est"]
cat(sprintf("    间接效应 = 总效应 - 直接效应 = %.3f  [真实 0.30]\n", indirect))
cat(sprintf("    中介比例 = %.1f%%  [真实 50%%]\n", prop * 100))

# 交叉验证: a*b 乘积法（a=X→M 用X工具; b=MVMR 校正后 M 效应）
a_est <- ivw(beta_x[1:Kx], beta_m[1:Kx], se_m[1:Kx])
ind_prod <- a_est["est"] * res_mv@Estimate[2]
cat(sprintf("    交叉验证 a*b 乘积法: a=%.3f, b=%.3f, 间接效应=%.3f\n",
            a_est["est"], res_mv@Estimate[2], ind_prod))

# ---------- 6. 保存 ----------
res <- data.frame(
  参数 = c("总效应", "直接效应(MVMR校正)", "间接效应(差值法)", "中介比例", "间接效应(乘积法)"),
  估计 = c(total["est"], direct, indirect, prop, ind_prod),
  真实值 = c(0.60, 0.30, 0.30, 0.50, 0.30))
print(res, row.names = FALSE)
write.csv(res, "02.analysis/mvmr_mediation/results.csv", row.names = FALSE)
cat("\n多变量中介 MR 完成 ✔\n")
