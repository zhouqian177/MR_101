#!/usr/bin/env Rscript
# 43_conmix_mr.R — Contamination mixture MR（ConMix 多效性稳健方法）
#
# ConMix（Bowden et al. 2016）对工具变量集建模为"有效"与"无效"两类，
# 有效 IV 效应为 θ，无效 IV 的效应为 θ + 多效性偏差。通过混合估计
# 同时识别有效 IV 子集与因果效应，对多效性高度稳健。
#
# 数据: LDL-C -> CHD 在线数据（28 个工具变量）
# 对比: IVW / Egger / 加权中位数 / ConMix
# 输出: 02.analysis/conmix/ 结果

options(width = 150)
dir.create("02.analysis/conmix", showWarnings = FALSE)
suppressMessages(library(MendelianRandomization))

cat("========================================\n")
cat("ConMix 多效性稳健方法\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 数据 ----------
cat("[1] 数据: LDL-C -> CHD（在线, 28 个工具变量）\n")
d <- read.csv("02.analysis/opengwas/online/harmonised.csv", stringsAsFactors = FALSE)
d <- subset(d, mr_keep)
cat("    工具变量数:", nrow(d), "\n")

mrin <- mr_input(bx = d$beta.exposure, bxse = d$se.exposure,
                 by = d$beta.outcome, byse = d$se.outcome,
                 snps = d$SNP, exposure = "LDL-C", outcome = "CHD")

# ---------- 各方法 ----------
cat("\n[2] 各方法对比\n")
# IVW
ivw <- mr_ivw(mrin)
cat(sprintf("    IVW:      beta=%.4f (SE=%.4f, P=%.3g)\n",
            ivw@Estimate, ivw@StdError, ivw@Pvalue))
# Egger
egg <- mr_egger(mrin)
cat(sprintf("    Egger:    beta=%.4f (SE=%.4f, P=%.3g), 截距 P=%.3g\n",
            egg@Estimate, egg@StdError.Est, egg@Pvalue.Est, egg@Pvalue.Int))
# 加权中位数
med <- mr_median(mrin)
cat(sprintf("    WMedian:  beta=%.4f (SE=%.4f, P=%.3g)\n",
            med@Estimate, med@StdError, med@Pvalue))
# ConMix
cat("\n[3] ConMix 分析（有效 IV 比例由混合模型自动估计）\n")
con <- mr_conmix(mrin)
cat(sprintf("    ConMix:   beta=%.4f (95%%CI: %.4f-%.4f, P=%.3g)\n",
            con@Estimate, con@CILower, con@CIUpper, con@Pvalue))
cat("    有效 IV 比例：基于 CIR 方法自动选择\n")

# ---------- 汇总对比 ----------
cat("\n[4] 结果汇总\n")
res <- data.frame(
  方法 = c("IVW", "Egger", "WMedian", "ConMix"),
  beta = c(ivw@Estimate, egg@Estimate, med@Estimate, con@Estimate),
  CI_lower = c(ivw@CILower, egg@CILower.Est, med@CILower, con@CILower),
  CI_upper = c(ivw@CIUpper, egg@CIUpper.Est, med@CIUpper, con@CIUpper),
  P = c(ivw@Pvalue, egg@Pvalue.Est, med@Pvalue, con@Pvalue),
  OR = exp(c(ivw@Estimate, egg@Estimate, med@Estimate, con@Estimate)))
print(res, row.names = FALSE)
write.csv(res, "02.analysis/conmix/results.csv", row.names = FALSE)

cat("\n解读: ConMix 通过混合模型区分有效/无效 IV，\n")
cat("      对多效性稳健。与 IVW/Egger/中位数对比可获得更全面的结论。\n")
cat("\nConMix 分析完成 ✔\n")