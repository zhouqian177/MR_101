#!/usr/bin/env Rscript
# 14_mrmix.R — MRMix 稳健混合模型（Qi & Chatterjee 2019）
#
# MRMix 假设工具变量分为"有效 IV"（无多效性）与"无效 IV"（存在水平多效性）
# 两类，用高斯混合模型建模，在存在水平多效性时仍可一致估计因果效应。
# 输出参数:
#   theta  - 因果效应估计
#   pi0    - 有效工具变量比例（pi0 越接近 1 多效性越弱）
#   sigma2 - 多效性方差
#
# 数据 1: 包内置 sumstats（388 SNP，直接可用）
# 数据 2: BMI15(暴露) -> MDD18(结局) 真实 GWAS 汇总数据应用
#
# 输出: 02.analysis/mrmix/ 下结果

options(width = 150)
dir.create("02.analysis/mrmix", showWarnings = FALSE)
suppressMessages(library(MRMix))

cat("========================================\n")
cat("MRMix 稳健混合模型\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 数据 1: 内置 sumstats ----------
cat("[1] 内置 sumstats 数据（388 SNP）\n")
data("sumstats", package = "MRMix")
r1 <- MRMix(sumstats$betahat_x, sumstats$betahat_y, sumstats$sx, sumstats$sy)
cat(sprintf("    theta(因果效应) = %.4f (SE=%.4f, P=%.3g)\n",
            r1$theta, r1$SE_theta, r1$pvalue_theta))
cat(sprintf("    pi0(有效IV比例) = %.3f | sigma2(多效性方差) = %.5f\n",
            r1$pi0, r1$sigma2))

# ---------- 数据 2: BMI15 -> MDD18 真实数据 ----------
cat("\n[2] 真实数据: BMI15(暴露) -> MDD18(结局)\n")
data("BMI15", package = "MRMix"); data("MDD18", package = "MRMix")
bmi <- BMI15[, c("SNP", "beta", "se")]
mdd <- MDD18[, c("SNP", "OR", "SE")]
mdd$beta_out <- log(mdd$OR)               # OR -> log OR
m <- merge(bmi, mdd, by = "SNP")
m <- m[complete.cases(m), ]
cat("    合并后 SNP 数:", nrow(m), "\n")
# 注: 教学演示未做等位基因方向对齐；完整流程需按效应等位基因 harmonise
r2 <- MRMix(m$beta, m$beta_out, m$se, m$SE)
cat(sprintf("    theta(BMI->MDD) = %.4f (SE=%.4f, P=%.3g)\n",
            r2$theta, r2$SE_theta, r2$pvalue_theta))
cat(sprintf("    pi0 = %.3f | sigma2 = %.5f\n", r2$pi0, r2$sigma2))

# ---------- 与 IVW 对比 ----------
cat("\n[3] 与 IVW 对比（内置数据）\n")
w <- 1 / sumstats$sy^2
ivw_b <- sum(sumstats$betahat_x * sumstats$betahat_y * w) / sum(sumstats$betahat_x^2 * w)
cat(sprintf("    IVW = %.4f | MRMix = %.4f\n", ivw_b, r1$theta))

res <- data.frame(
  数据集 = c("内置sumstats", "BMI15->MDD18"),
  theta = c(r1$theta, r2$theta),
  SE = c(r1$SE_theta, r2$SE_theta),
  P = c(r1$pvalue_theta, r2$pvalue_theta),
  pi0 = c(r1$pi0, r2$pi0),
  sigma2 = c(r1$sigma2, r2$sigma2))
write.csv(res, "02.analysis/mrmix/results.csv", row.names = FALSE)
print(res, row.names = FALSE)
cat("\nMRMix 完成 ✔\n")
