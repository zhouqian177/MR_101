#!/usr/bin/env Rscript
# 00_env.R — MR 分析环境检查：R 版本、关键 R 包、系统工具
# 用法: Rscript scripts/00_env.R
# 输出: 01.tools/env_check.txt（同时打印到终端）

options(width = 120)
out <- file("01.tools/env_check.txt", open = "w")
sink(out, type = "output")

cat("========================================\n")
cat("MR_101 环境检查报告\n")
cat("检查时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

cat("[1] R 版本\n")
cat("    ", R.version.string, "\n")
cat("    平台:", R.version$platform, "\n\n")

cat("[2] 关键 R 包\n")
pkgs <- c("TwoSampleMR", "MendelianRandomization", "ieugwasr",
          "MRInstruments", "MRPRESSO", "data.table", "ggplot2",
          "dplyr", "knitr", "rmarkdown")
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  if (ok) {
    v <- as.character(packageVersion(p))
    cat(sprintf("    %-22s OK     v%s\n", p, v))
  } else {
    cat(sprintf("    %-22s MISSING\n", p))
  }
}

cat("\n[3] 系统工具\n")
sys <- c("Rscript", "python3", "plink", "plink2", "gcta", "bcftools", "samtools")
for (t in sys) {
  w <- Sys.which(t)
  if (nzchar(w)) cat(sprintf("    %-10s %s\n", t, w)) else cat(sprintf("    %-10s (未找到，非必需)\n", t))
}

cat("\n[4] OpenGWAS API 连通性\n")
net <- tryCatch({
  suppressMessages(requireNamespace("ieugwasr", quietly = TRUE))
  if (requireNamespace("ieugwasr", quietly = TRUE)) {
    x <- ieugwasr::gwasinfo("ieu-b-110")
    sprintf("    连通 OK，ieu-b-110 样本量=%s\n", x$ncase[1])
  } else "    ieugwasr 未安装，跳过\n"
}, error = function(e) sprintf("    连接失败: %s\n", conditionMessage(e)))
cat(net, "\n")

cat("========================================\n")
cat("环境检查完成\n")
sink()
close(out)
cat("报告已写入: 01.tools/env_check.txt\n")
