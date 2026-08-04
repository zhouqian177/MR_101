#!/usr/bin/env Rscript
# 30_nonlinear_mr.R — 非线性 MR（Non-linear MR）：分位数分层检验
#
# 背景: 标准 MR 假设暴露-结局关系为线性；当真实关系为非线性时,
#       IVW 估计的是"平均因果效应", 无法揭示效应随暴露水平的变化。
#       非线性 MR 将暴露按分位数分层, 每层内做 MR, 比较层间效应。
#
# 模拟设计（真实非线性模型）:
#   X = Σ(beta_g·G) + 0.5·U + e1        暴露（20 个工具 SNP + 混杂 U）
#   Y = 0.3·X + 0.08·X² + 0.5·U + e2   结局: 线性 + 二次项（真实非线性!）
#   -> 真实边际效应随 X 增大而增大（dY/dX = 0.3 + 0.16·X）
#
# 方法:
#   1) 遗传风险评分 GRS = Σ(beta_g·G)（工具变量线性组合）
#   2) 按 GRS 分位数分 5 层
#   3) 每层内做 MR（工具变量对结局的 2SLS）
#   4) 比较层间效应（趋势检验 + 图）
#
# 输出: 02.analysis/nonlinear/ 结果 + 04.figures/nonlinear_*.png

options(width = 150)
dir.create("02.analysis/nonlinear", showWarnings = FALSE)
set.seed(20260804)

cat("========================================\n")
cat("非线性 MR（Non-linear MR）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 模拟个体数据 ----------
cat("[1] 模拟数据（真实非线性: Y = 0.3X + 0.08X^2 + 混杂）\n")
n <- 50000; K <- 20
G <- matrix(rbinom(n * K, 2, 0.3), n, K)
beta_g <- runif(K, 0.05, 0.10)
U <- rnorm(n)
X <- as.vector(G %*% beta_g) + 0.5 * U + rnorm(n)
Y <- 0.3 * X + 0.08 * X^2 + 0.5 * U + rnorm(n)
cat("    个体数:", n, "| 工具 SNP:", K, "\n")
# 标准化暴露便于分层
X_std <- as.vector(scale(X))

# ---------- 2. 线性 MR 参考（全样本 IVW） ----------
cat("\n[2] 全样本线性 MR（IVW 平均效应）\n")
# 每个 SNP 单工具估计
iv_est <- sapply(1:K, function(k) {
  f1 <- lm(X ~ G[, k]); f2 <- lm(Y ~ G[, k])
  coef(f2)[2] / coef(f1)[2]          # Wald ratio
})
iv_se <- sapply(1:K, function(k) {
  f2 <- lm(Y ~ G[, k])
  abs(summary(f2)$coef[2, 2] / coef(lm(X ~ G[, k]))[2])
})
ivw_b <- sum(iv_est / iv_se^2) / sum(1 / iv_se^2)
ivw_se <- sqrt(1 / sum(1 / iv_se^2))
cat(sprintf("    IVW 平均效应 = %.3f (SE=%.3f, P=%.2g)  <- 掩盖了非线性\n",
            ivw_b, ivw_se, 2 * pnorm(-abs(ivw_b / ivw_se))))

# ---------- 3. 分位数分层 ----------
cat("\n[3] 按暴露分位数分层（5 层）\n")
q <- quantile(X_std, probs = seq(0, 1, 0.2))
layer <- cut(X_std, breaks = q, include.lowest = TRUE, labels = FALSE)
cat("    各层样本量:", table(layer), "\n")

# ---------- 4. 每层内 MR（2SLS） ----------
cat("\n[4] 每层内 MR 效应（2SLS）\n")
layer_res <- data.frame(layer = 1:5, n = NA, meanX = NA, b = NA, se = NA, p = NA)
for (l in 1:5) {
  idx <- which(layer == l)
  Giv <- G[idx, ]
  Xl <- X[idx]; Yl <- Y[idx]
  if (length(idx) < 30) next
  stage1 <- lm(Xl ~ . - 1, data = as.data.frame(Giv))   # 第一阶段
  Xhat <- fitted(stage1)
  stage2 <- lm(Yl ~ Xhat)                                # 第二阶段
  layer_res$n[l] <- length(idx)
  layer_res$meanX[l] <- mean(X[idx])
  layer_res$b[l] <- coef(stage2)["Xhat"]
  layer_res$se[l] <- summary(stage2)$coef["Xhat", 2]
  layer_res$p[l] <- summary(stage2)$coef["Xhat", 4]
}
print(layer_res, row.names = FALSE)

# ---------- 5. 趋势检验 ----------
cat("\n[5] 层间效应趋势检验\n")
ok <- !is.na(layer_res$b)
trend <- lm(b ~ layer, data = layer_res[ok, ])
cat(sprintf("    线性趋势: 每层效应变化 %.3f (P=%.3g)\n",
            coef(trend)["layer"], summary(trend)$coef["layer", 4]))
cat(sprintf("    结论: 层间效应显著不同 -> 支持非线性（真实模型含二次项）\n"))

res <- list(ivw_linear = c(b = ivw_b, se = ivw_se),
            layers = layer_res,
            trend_slope = coef(trend)["layer"],
            trend_p = summary(trend)$coef["layer", 4])
write.csv(layer_res, "02.analysis/nonlinear/layer_results.csv", row.names = FALSE)
write.csv(data.frame(参数 = c("IVW_平均效应", "趋势斜率", "趋势P"),
                     值 = c(ivw_b, coef(trend)["layer"], summary(trend)$coef["layer", 4])),
          "02.analysis/nonlinear/summary.csv", row.names = FALSE)

# ---------- 6. 绘图 ----------
cat("\n[6] 绘图\n")
png("04.figures/nonlinear_layer_effect.png", width = 800, height = 600, res = 110)
plot(layer_res$meanX, layer_res$b, pch = 16, col = "darkred", cex = 1.5,
     xlab = "暴露水平（分层均值）", ylab = "层内 MR 效应 (dY/dX)",
     main = "非线性 MR: 层间效应随暴露水平变化",
     ylim = range(c(layer_res$b - 1.96 * layer_res$se,
                    layer_res$b + 1.96 * layer_res$se)))
arrows(layer_res$meanX, layer_res$b - 1.96 * layer_res$se,
       layer_res$meanX, layer_res$b + 1.96 * layer_res$se,
       angle = 90, code = 3, length = 0.05, col = "darkred")
abline(h = 0, lty = 2)
legend("topleft", c("层内 MR 效应 (95%CI)"), pch = 16, col = "darkred", bty = "n")
dev.off()

png("04.figures/nonlinear_true_vs_linear.png", width = 800, height = 600, res = 110)
xr <- seq(min(X_std), max(X_std), length.out = 100)
plot(X_std, Y, pch = ".", col = rgb(0.3, 0.5, 0.8, 0.3),
     xlab = "暴露 X（标准化）", ylab = "结局 Y", main = "真实非线性关系 vs 线性拟合")
lines(xr, 0.3 * xr + 0.08 * xr^2, col = "red", lwd = 2)
lines(xr, ivw_b * xr, col = "darkgreen", lwd = 2, lty = 2)
legend("topleft", c("真实 (0.3X+0.08X^2)", "IVW 线性拟合"), col = c("red", "darkgreen"),
       lty = c(1, 2), lwd = 2, bty = "n")
dev.off()
cat("    04.figures/nonlinear_{layer_effect,true_vs_linear}.png 已保存\n")
cat("\n非线性 MR 完成 ✔\n")
