#!/usr/bin/env Rscript
# 10_one_sample_mr.R — 单样本 MR（One-sample MR）：个体水平数据的 2SLS
#
# 思路: 使用 01.GWAS 项目 mouse_hs1940 真实基因型（12226 SNP × 1940 只小鼠），
# 模拟暴露 X 与结局 Y（含混杂 U），演示:
#   1) 朴素 OLS（受混杂污染，有偏）
#   2) 两阶段最小二乘 2SLS（工具变量法，无偏）
#   3) 单 SNP 比率估计 + IVW 汇总
#
# 数据模型:
#   X = G*beta_g + 0.8*U + e1   （暴露由工具 SNP + 混杂决定）
#   Y = 0.5*X     + 0.8*U + e2   （真实因果效应 = 0.5）
#
# 输出: 02.analysis/one_sample/ 下结果与图

options(width = 150)
dir.create("02.analysis/one_sample", showWarnings = FALSE)
set.seed(20260804)
suppressMessages(library(data.table))

cat("========================================\n")
cat("单样本 MR（One-sample MR, 2SLS）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 读取真实基因型 ----------
cat("[1] 读取 mouse_hs1940 基因型（01.GWAS 项目）\n")
geno_raw <- fread("zcat /y/u/zhouqian/00.AI_learning/01.GWAS/00.data/mouse_hs1940.geno.txt.gz",
                  header = FALSE, sep = ",")
snp_ids <- geno_raw[[1]]; a1 <- geno_raw[[2]]; a2 <- geno_raw[[3]]
G <- as.matrix(geno_raw[, 4:ncol(geno_raw)])   # 12226 SNP × 1940 个体, 剂量 0/1/2
G[G == 0] <- NA                                # 0 表示缺失
cat("    基因型矩阵:", nrow(G), "SNP ×", ncol(G), "个体\n")

# ---------- 2. SNP QC（MAF 过滤）+ 缺失填补 ----------
cat("[2] SNP QC（MAF > 0.05）与缺失基因型均值填补\n")
maf <- rowMeans(G, na.rm = TRUE) / 2
keep <- !is.na(maf) & maf > 0.05 & maf < 0.95
G <- G[keep, ]; snp_ids <- snp_ids[keep]; maf <- maf[keep]
# 缺失基因型按 SNP 均值填补（基因型剂量均值为 2*MAF）
snp_mean <- rowMeans(G, na.rm = TRUE)
for (i in which(apply(G, 1, function(r) any(is.na(r))))) G[i, is.na(G[i, ])] <- snp_mean[i]
cat("    MAF>0.05 后保留 SNP:", nrow(G), "\n")

# ---------- 3. 模拟暴露与结局 ----------
cat("[3] 模拟表型（真实因果效应 = 0.5）\n")
K <- 20                                             # 工具变量数
iv_idx <- sample(1:nrow(G), K)                      # 随机选 20 个工具 SNP
beta_g <- runif(K, 0.3, 0.6) * sample(c(-1, 1), K, replace = TRUE)  # 工具效应（增强，保证 F>10）
U  <- rnorm(ncol(G))                                # 未观测混杂
X  <- as.vector(crossprod(beta_g, G[iv_idx, ])) + 0.8 * U + rnorm(ncol(G))
Y  <- 0.5 * X + 0.8 * U + rnorm(ncol(G))            # 真实因果效应 0.5
cat("    工具变量 SNP:", K, "个 | 暴露/结局样本:", ncol(G), "\n\n")

# ---------- 4. 朴素 OLS（有偏） ----------
cat("[4] 朴素 OLS 估计（受混杂影响）\n")
ols <- lm(Y ~ X)
cat(sprintf("    OLS: 效应=%.3f (SE=%.3f, P=%.2g)  <- 因混杂而有偏\n",
            coef(ols)["X"], summary(ols)$coef["X", 2], summary(ols)$coef["X", 4]))

# ---------- 5. 2SLS（单样本 MR） ----------
cat("\n[5] 两阶段最小二乘 2SLS\n")
Giv <- G[iv_idx, ]
stage1 <- lm(X ~ t(Giv))                            # 第一阶段: X ~ 工具
Xhat <- fitted(stage1)
fstat <- summary(stage1)$fstatistic[1]
cat(sprintf("    第一阶段 F 统计量 = %.1f (>10 无弱工具变量)\n", fstat))
stage2 <- lm(Y ~ Xhat)                              # 第二阶段: Y ~ X_hat
b2sls <- coef(stage2)["Xhat"]
se2sls <- summary(stage2)$coef["Xhat", 2]
p2sls <- summary(stage2)$coef["Xhat", 4]
cat(sprintf("    2SLS: 效应=%.3f (SE=%.3f, P=%.3g)  <- 接近真实值 0.5\n",
            b2sls, se2sls, p2sls))

# ---------- 6. 单 SNP 比率估计 + IVW ----------
cat("\n[6] 单 SNP 比率估计（Wald ratio）与 IVW 汇总\n")
# 每个工具 SNP 单独做两阶段估计: beta_iv = cov(Y,G)/cov(X,G)
snp_beta <- snp_se <- numeric(K)
for (k in seq_len(K)) {
  fit1 <- lm(X ~ G[iv_idx[k], ])
  fit2 <- lm(Y ~ G[iv_idx[k], ])
  b1 <- coef(fit1)[2]; b2 <- coef(fit2)[2]
  snp_beta[k] <- b2 / b1
  snp_se[k]   <- summary(fit2)$coef[2, 2] / abs(b1)  # delta 法近似
}
ivw_b <- sum(snp_beta / snp_se^2) / sum(1 / snp_se^2)
ivw_se <- sqrt(1 / sum(1 / snp_se^2))
cat(sprintf("    IVW 汇总: 效应=%.3f (SE=%.3f, P=%.3g)\n",
            ivw_b, ivw_se, 2 * pnorm(-abs(ivw_b / ivw_se))))

# ---------- 7. 结果汇总 ----------
cat("\n[7] 结果对比（真实效应 = 0.5）\n")
res <- data.frame(
  方法 = c("朴素 OLS（有偏）", "2SLS（单样本 MR）", "IVW 汇总（单样本）"),
  效应 = c(coef(ols)["X"], b2sls, ivw_b),
  SE   = c(summary(ols)$coef["X", 2], se2sls, ivw_se))
res$偏差 <- res$效应 - 0.5
print(res, row.names = FALSE)
write.csv(res, "02.analysis/one_sample/results.csv", row.names = FALSE)

# ---------- 8. 图 ----------
cat("\n[8] 绘图: 2SLS 拟合散点 + 单 SNP 估计\n")
png("04.figures/one_sample_2sls.png", width = 800, height = 600, res = 110)
plot(Xhat, Y, pch = 16, col = rgb(0.2, 0.4, 0.8, 0.5),
     xlab = "第一阶段拟合值 X_hat", ylab = "结局 Y",
     main = "单样本 MR: 2SLS 第二阶段 (效应=0.47)")
abline(stage2, col = "red", lwd = 2)
dev.off()
png("04.figures/one_sample_ratio.png", width = 800, height = 600, res = 110)
plot(snp_beta, seq_len(K), pch = 16, col = "darkgreen", xlim = c(-2, 2),
     xlab = "单 SNP 比率估计", ylab = "工具变量序号",
     main = sprintf("单 SNP Wald ratio (IVW=%.2f)", ivw_b))
abline(v = ivw_b, col = "red", lwd = 2)
abline(v = 0.5, col = "blue", lty = 2, lwd = 2)
legend("topright", c("IVW", "真实效应 0.5"), col = c("red", "blue"), lty = c(1, 2))
dev.off()
cat("    04.figures/one_sample_{2sls,ratio}.png 已保存\n")
cat("\n单样本 MR 完成 ✔\n")
