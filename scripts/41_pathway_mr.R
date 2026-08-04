#!/usr/bin/env Rscript
# 41_pathway_mr.R — 通路 MR（Pathway MR / Sequential Mediation）
#
# 背景: 两步中介（扩展3）只考虑一个中介; 通路 MR 扩展为多步中介链:
#   X -> M1 -> M2 -> Y
#   - 间接效应(路径1): X -> M1 -> Y（单中介）
#   - 间接效应(路径2): X -> M1 -> M2 -> Y（双中介串联）
#   - 总间接效应 = 路径1 + 路径2
#   - 直接效应 = 总效应 - 总间接效应
#
# 模拟设计:
#   X = Gx*bx + e1               暴露（15 个工具 SNP）
#   M1 = a1*X + G1*b1 + e2       中介1（X→M1, a1=0.5）
#   M2 = a2*M1 + G2*b2 + e3      中介2（M1→M2, a2=0.6）
#   Y = c'*X + b1*M1 + b2*M2 + e4 结局（c'=0.2, b1=0.3, b2=0.4）
#
# 真实值:
#   总效应 = 0.2 + 0.5*0.3 + 0.5*0.6*0.4 = 0.2+0.15+0.12 = 0.47
#   直接效应 = 0.20
#   间接(单中介) = 0.5*0.3 = 0.15
#   间接(双中介) = 0.5*0.6*0.4 = 0.12
#   总间接 = 0.27；中介比例 = 0.27/0.47 = 57.4%
#
# 输出: 02.analysis/pathway/ 结果

options(width = 150)
dir.create("02.analysis/pathway", showWarnings = FALSE)
set.seed(20260804)

cat("========================================\n")
cat("通路 MR（Pathway MR / Sequential Mediation）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 模拟个体数据 ----------
cat("[1] 模拟数据（X→M1→M2→Y）\n")
n <- 50000; Kx <- 15; K1 <- 15; K2 <- 15
Gx <- matrix(rbinom(n * Kx, 2, 0.3), n, Kx)
G1 <- matrix(rbinom(n * K1, 2, 0.3), n, K1)
G2 <- matrix(rbinom(n * K2, 2, 0.3), n, K2)
bx <- runif(Kx, 0.06, 0.12); b1 <- runif(K1, 0.06, 0.12); b2 <- runif(K2, 0.06, 0.12)
X <- as.vector(Gx %*% bx) + rnorm(n)
M1 <- 0.5 * X + as.vector(G1 %*% b1) + rnorm(n)
M2 <- 0.6 * M1 + as.vector(G2 %*% b2) + rnorm(n)
Y <- 0.2 * X + 0.3 * M1 + 0.4 * M2 + rnorm(n)
cat("    个体数:", n, "| 工具: X", Kx, "SNP, M1", K1, "SNP, M2", K2, "SNP\n")

# ---------- 2. 生成汇总统计量 ----------
cat("[2] 生成汇总统计量\n")
gen_sumstats <- function(G, Y) {
  K <- ncol(G); b <- se <- numeric(K)
  for (k in 1:K) { f <- summary(lm(Y ~ G[, k])); b[k] <- f$coef[2, 1]; se[k] <- f$coef[2, 2] }
  list(b = b, se = se)
}
bx <- gen_sumstats(Gx, X); sex <- bx$se
bx_m1 <- gen_sumstats(Gx, M1)
b1_m2 <- gen_sumstats(G1, M2)
b2_y <- gen_sumstats(G2, Y)
bx_y <- gen_sumstats(Gx, Y)

# ---------- IVW 函数 ----------
ivw <- function(bx, by, sey) {
  w <- 1 / sey^2; b <- sum(bx * by * w) / sum(bx^2 * w)
  se <- sqrt(1 / sum(bx^2 * w)); c(est = b, se = se, p = 2 * pnorm(-abs(b / se)))
}

# ---------- 3. 路径分解 ----------
cat("\n[3] 路径分解\n")
# 总效应: X -> Y
total <- ivw(bx$b, bx_y$b, bx_y$se)
cat(sprintf("    总效应 X->Y          = %.3f (SE=%.3f)  [真实 0.47]\n", total["est"], total["se"]))
# 路径1: X -> M1 -> Y
a1 <- ivw(bx$b, bx_m1$b, bx_m1$se)
b1_m2y <- ivw(b1_m2$b, b2_y$b, b2_y$se)  # M1 -> Y (用 M1 工具)
# 路径2: X -> M1 -> M2 -> Y
a2 <- ivw(gen_sumstats(G1, M1)$b, gen_sumstats(G1, M2)$b, gen_sumstats(G1, M2)$se)  # M1->M2 (用M1工具)
b2 <- ivw(b2_y$b, b2_y$b, b2_y$se)  # 注意: M2->Y 用 M2 工具

# 更准确的: 需要单独的 M2->Y 估计
# 直接用 M2 工具对 Y 回归
b2_m2y <- gen_sumstats(G2, M2)$b  # M2 工具对 M2 的效应
b2_b_y <- gen_sumstats(G2, Y)$b    # M2 工具对 Y 的效应
se2_b_y <- gen_sumstats(G2, Y)$se
path_M2 <- ivw(b2_m2y, b2_b_y, se2_b_y)
cat(sprintf("    a1 X->M1             = %.3f (SE=%.3f)  [真实 0.50]\n", a1["est"], a1["se"]))
cat(sprintf("    path M1->Y(单中介)   = %.3f (SE=%.3f)  [真实 0.30]\n", b1_m2y["est"], b1_m2y["se"]))
cat(sprintf("    path M2->Y(校正后)   = %.3f (SE=%.3f)  [真实 0.40]\n", path_M2["est"], path_M2["se"]))

# 间接效应分解
ind_single <- a1["est"] * b1_m2y["est"]  # X->M1->Y
ind_double <- a1["est"] * a2["est"] * path_M2["est"]  # X->M1->M2->Y
ind_total <- ind_single + ind_double
direct <- total["est"] - ind_total
prop <- ind_total / total["est"]
cat(sprintf("\n    间接效应(单中介 X->M1->Y) = %.3f  [真实 0.15]\n", ind_single))
cat(sprintf("    间接效应(双中介串联)     = %.3f  [真实 0.12]\n", ind_double))
cat(sprintf("    总间接效应               = %.3f  [真实 0.27]\n", ind_total))
cat(sprintf("    直接效应                 = %.3f  [真实 0.20]\n", direct))
cat(sprintf("    中介比例                 = %.1f%% [真实 57.4%%]\n", prop * 100))

res <- data.frame(
  参数 = c("总效应", "a1(X->M1)", "单中介(M1->Y)", "双中介(M2->Y)",
           "间接(单中介)", "间接(双中介)", "总间接", "直接", "中介比例"),
  估计 = c(total["est"], a1["est"], b1_m2y["est"], path_M2["est"],
           ind_single, ind_double, ind_total, direct, prop),
  真实值 = c(0.47, 0.50, 0.30, 0.40, 0.15, 0.12, 0.27, 0.20, 0.574))
write.csv(res, "02.analysis/pathway/results.csv", row.names = FALSE)
print(res, row.names = FALSE)
cat("\n通路 MR 完成 ✔\n")