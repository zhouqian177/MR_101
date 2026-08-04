#!/usr/bin/env Rscript
# 12_mediation_mr.R — 两步中介 MR（Two-step Mediation MR）
#
# 模型（真实参数）:
#   X = Gx*bx + e1                   暴露（15 个工具 SNP）
#   M = 0.5*X + Gm*bm + e2           中介（受暴露影响, a=0.5）
#   Y = 0.3*X + 0.6*M + e3           结局（直接效应 c'=0.3, b=0.6）
#
# 推导: 间接效应 = a*b = 0.30；总效应 = c' + a*b = 0.60；中介比例 = 0.50
#
# 两步 MR:
#   第一步: X -> M  的因果效应 a（用 X 的工具 SNP, IVW）
#   第二步: M -> Y  的因果效应 b（用 M 的工具 SNP, IVW）
#   间接效应 = a*b；直接效应 = 总效应 - 间接效应
#
# 输出: 02.analysis/mediation/ 下结果

options(width = 150)
dir.create("02.analysis/mediation", showWarnings = FALSE)
set.seed(20260804)

cat("========================================\n")
cat("两步中介 MR（Two-step Mediation MR）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 模拟个体数据 ----------
cat("[1] 模拟数据（a=0.5, b=0.6, c'=0.3）\n")
n <- 30000
Kx <- 15; Km <- 15
Gx <- matrix(rbinom(n * Kx, 2, 0.3), n, Kx)   # 暴露工具
Gm <- matrix(rbinom(n * Km, 2, 0.3), n, Km)   # 中介工具
bx <- runif(Kx, 0.06, 0.12); bm <- runif(Km, 0.06, 0.12)
X <- as.vector(Gx %*% bx) + rnorm(n)
M <- 0.5 * X + as.vector(Gm %*% bm) + rnorm(n)
Y <- 0.3 * X + 0.6 * M + rnorm(n)
cat("    个体数:", n, "| 工具: X 用", Kx, "SNP, M 用", Km, "SNP\n")

# ---------- 2. 生成汇总统计量 ----------
cat("[2] 生成汇总统计量（每个 SNP 对 X/M/Y 的 beta, se）\n")
beta_x <- beta_m <- beta_y <- se_x <- se_m <- se_y <- numeric(Kx + Km)
for (k in 1:Kx) {
  g <- Gx[, k]
  beta_x[k] <- coef(summary(lm(X ~ g)))[2, 1]; se_x[k] <- coef(summary(lm(X ~ g)))[2, 2]
  beta_m[k] <- coef(summary(lm(M ~ g)))[2, 1]; se_m[k] <- coef(summary(lm(M ~ g)))[2, 2]
  beta_y[k] <- coef(summary(lm(Y ~ g)))[2, 1]; se_y[k] <- coef(summary(lm(Y ~ g)))[2, 2]
}
for (k in 1:Km) {
  g <- Gm[, k]
  j <- Kx + k
  beta_x[j] <- coef(summary(lm(X ~ g)))[2, 1]; se_x[j] <- coef(summary(lm(X ~ g)))[2, 2]
  beta_m[j] <- coef(summary(lm(M ~ g)))[2, 1]; se_m[j] <- coef(summary(lm(M ~ g)))[2, 2]
  beta_y[j] <- coef(summary(lm(Y ~ g)))[2, 1]; se_y[j] <- coef(summary(lm(Y ~ g)))[2, 2]
}
cat("    汇总统计量:", Kx + Km, "SNP\n")

# ---------- 3. IVW 辅助函数 ----------
ivw <- function(bx, by, byse) {
  w <- 1 / byse^2
  b <- sum(bx * by * w) / sum(bx^2 * w)
  se <- sqrt(1 / sum(bx^2 * w))
  c(est = b, se = se, p = 2 * pnorm(-abs(b / se)))
}

# ---------- 4. 两步 MR ----------
cat("\n[3] 两步 MR\n")
# 第一步: X -> M（用 X 的工具 SNP）
step1 <- ivw(beta_x[1:Kx], beta_m[1:Kx], se_m[1:Kx])
cat(sprintf("    第一步 X->M : a = %.3f (SE=%.3f, P=%.2g)  [真实 0.5]\n",
            step1["est"], step1["se"], step1["p"]))
# 第二步: M -> Y（用 M 的工具 SNP）
step2 <- ivw(beta_m[Kx + 1:Km], beta_y[Kx + 1:Km], se_y[Kx + 1:Km])
cat(sprintf("    第二步 M->Y : b = %.3f (SE=%.3f, P=%.2g)  [真实 0.6]\n",
            step2["est"], step2["se"], step2["p"]))
# 总效应: X -> Y（用 X 的工具 SNP）
total <- ivw(beta_x[1:Kx], beta_y[1:Kx], se_y[1:Kx])
cat(sprintf("    总效应 X->Y : %.3f (SE=%.3f)  [真实 0.60]\n", total["est"], total["se"]))

# ---------- 5. 中介效应分解 ----------
cat("\n[4] 中介效应分解\n")
indirect <- step1["est"] * step2["est"]
direct   <- total["est"] - indirect
prop     <- indirect / total["est"]
# delta 法近似间接效应 SE: se(a*b) ≈ sqrt(b^2*se_a^2 + a^2*se_b^2)
se_ind  <- sqrt(step2["est"]^2 * step1["se"]^2 + step1["est"]^2 * step2["se"]^2)
cat(sprintf("    间接效应 a*b     = %.3f (SE=%.3f)  [真实 0.30]\n", indirect, se_ind))
cat(sprintf("    直接效应 c'       = %.3f            [真实 0.30]\n", direct))
cat(sprintf("    中介比例         = %.1f%%           [真实 50%%]\n", prop * 100))

res <- data.frame(
  参数 = c("a (X->M)", "b (M->Y)", "总效应", "间接效应", "直接效应", "中介比例"),
  估计 = c(step1["est"], step2["est"], total["est"], indirect, direct, prop),
  真实值 = c(0.5, 0.6, 0.6, 0.3, 0.3, 0.5))
print(res, row.names = FALSE)
write.csv(res, "02.analysis/mediation/results.csv", row.names = FALSE)
cat("\n两步中介 MR 完成 ✔\n")
