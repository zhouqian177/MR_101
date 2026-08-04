#!/usr/bin/env Rscript
# 32_bayesian_mr.R — 贝叶斯 MR: cML（约束最大似然, Xue et al. 2021）
#
# 背景: cML 用约束最大似然方法, 允许工具变量存在相关/不相关的水平多效性,
#       通过 BIC/AIC 选择"有效工具变量"集合, 对多效性稳健。
#       cML-MA-BIC: 对所有候选模型做模型平均（更稳健）
#       cML-BIC-DP: 数据扰动版本, 校正选择不确定性
#
# 数据: LDL-C -> CHD 在线 harmonised 数据（28 个工具变量, 02.analysis/opengwas/online/harmonised.csv）
# 输出: 02.analysis/bayesian/ 结果 + 与 IVW/Egger 对比

options(width = 150)
dir.create("02.analysis/bayesian", showWarnings = FALSE)
suppressMessages(library(MendelianRandomization))

cat("========================================\n")
cat("贝叶斯 MR（cML 约束最大似然）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 数据 ----------
cat("[1] 数据: LDL-C -> CHD（OpenGWAS 在线, 28 个工具变量）\n")
d <- read.csv("02.analysis/opengwas/online/harmonised.csv", stringsAsFactors = FALSE)
d <- subset(d, mr_keep)
cat("    工具变量数:", nrow(d), "\n")

mrin <- mr_input(bx = d$beta.exposure, bxse = d$se.exposure,
                 by = d$beta.outcome, byse = d$se.outcome,
                 snps = d$SNP, exposure = "LDL-C", outcome = "CHD")

# ---------- 2. 参考方法 ----------
cat("\n[2] 参考方法（IVW / MR-Egger）\n")
ivw <- mr_ivw(mrin)
egger <- mr_egger(mrin)
cat(sprintf("    IVW : beta=%.4f (SE=%.4f, P=%.3g)\n",
            ivw@Estimate, ivw@StdError, ivw@Pvalue))
cat(sprintf("    Egger: beta=%.4f (SE=%.4f, P=%.3g), 截距 P=%.3g\n",
            egger@Estimate, egger@StdError.Est, egger@Pvalue.Est,
            egger@Pvalue.Int))

# ---------- 3. cML（约束最大似然） ----------
cat("\n[3] cML 分析（BIC 模型选择 + 模型平均 MA）\n")
cml <- tryCatch(mr_cML(mrin, n = 188578, random_seed = 20260804),
                error = function(e) paste("cML 失败:", conditionMessage(e)))
if (is.character(cml)) {
  cat("    ", cml, "\n")
} else {
  cat(sprintf("    cML-BIC(MA): beta=%.4f (SE=%.4f, P=%.3g)\n",
              cml@Estimate, cml@StdError, cml@Pvalue))
  cat("    cML 检出无效工具变量(BIC):", ifelse(length(cml@BIC_invalid) == 0, "无",
      paste(cml@BIC_invalid, collapse = ",")), "\n")
}

# ---------- 4. 汇总对比 ----------
cat("\n[4] 结果对比\n")
res <- data.frame(
  方法 = c("IVW", "MR-Egger", "cML"),
  beta = c(ivw@Estimate, egger@Estimate,
           ifelse(is.character(cml), NA, cml@Estimate)),
  SE = c(ivw@StdError, egger@StdError.Est,
         ifelse(is.character(cml), NA, cml@StdError)),
  P = c(ivw@Pvalue, egger@Pvalue.Est,
        ifelse(is.character(cml), NA, cml@Pvalue)))
res$OR <- exp(res$beta)
print(res, row.names = FALSE)
write.csv(res, "02.analysis/bayesian/results.csv", row.names = FALSE)

cat("\n解读: cML 在存在多效性时仍可稳健估计因果效应;\n")
cat("     与 IVW/Egger 对比可检验结论是否受多效性影响。\n")
cat("\n贝叶斯 MR（cML）完成 ✔\n")
