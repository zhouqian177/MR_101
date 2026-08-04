#!/usr/bin/env Rscript
# 03_sensitivity.R — MR 敏感性分析：MR-PRESSO、leave-one-out、单 SNP
# 输入: 02.analysis/harmonised.csv
# 输出: 02.analysis/sensitivity/ 下的 CSV 与 04.figures/ 图

options(width = 150)
dir.create("02.analysis/sensitivity", showWarnings = FALSE)
suppressMessages(library(TwoSampleMR))
suppressMessages(library(MRPRESSO))

cat("========================================\n")
cat("MR 敏感性分析（Telomere_length -> CHD）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

dat <- read.csv("02.analysis/harmonised.csv", stringsAsFactors = FALSE)
dat <- subset(dat, mr_keep)

# ---------- 1. 单 SNP 分析 ----------
cat("[1] 单 SNP 分析 + leave-one-out\n")
ss <- mr_singlesnp(dat)
loo <- mr_leaveoneout(dat)
write.csv(ss, "02.analysis/sensitivity/single_snp.csv", row.names = FALSE)
write.csv(loo, "02.analysis/sensitivity/leave_one_out.csv", row.names = FALSE)
cat("    单 SNP:", nrow(ss), "条 | leave-one-out:", nrow(loo), "条\n\n")

# ---------- 2. MR-PRESSO（离群值与水平多效性） ----------
cat("[2] MR-PRESSO（全局检验 + 离群值检测）\n")
presso <- tryCatch({
  d <- data.frame(BX = dat$beta.exposure, BY = dat$beta.outcome,
                  BXse = dat$se.exposure, BYse = dat$se.outcome)
  mr_presso(BetaOutcome = "BY", BetaExposure = "BX",
            SdOutcome = "BYse", SdExposure = "BXse", OUTLIERtest = TRUE,
            DISTORTIONtest = TRUE, data = d, NbDistribution = 1000)
}, error = function(e) paste("MR-PRESSO 失败:", conditionMessage(e)))

if (is.character(presso)) {
  cat("    ", presso, "\n")
} else {
  cat("    全局检验 P =", presso$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
  cat("    离群 SNP:", ifelse(length(presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`) == 0,
                              "无", paste(presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, collapse = ", ")), "\n")
  sink("02.analysis/sensitivity/mr_presso.txt")
  print(presso)
  sink()
  cat("    详情见 02.analysis/sensitivity/mr_presso.txt\n\n")
}

# ---------- 3. leave-one-out 图 ----------
cat("[3] leave-one-out 森林图\n")
ploo <- mr_leaveoneout_plot(loo)
ggplot2::ggsave("04.figures/mr_leaveoneout.png", ploo[[1]], width = 7, height = 8)
cat("    已保存 04.figures/mr_leaveoneout.png\n")

cat("\n敏感性分析完成 ✔\n")
