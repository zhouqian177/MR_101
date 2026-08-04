#!/usr/bin/env Rscript
# 20_opengwas_harmonise.R — 公开 GWAS 数据: LDL-C(暴露) 与 CHD(结局) harmonise
#
# 数据:
#   暴露: GLGC 2013 LDL (ieu-b-110 同源) -> 02.analysis/ldl_instruments_raw.txt (68 个工具变量)
#   结局: CARDIoGRAMplusC4D 2015 CAD (ieu-a-7 同源) -> 00.data/gwas/CHD_CAD2015.h.tsv.gz
#
# 流程: 按 rsid 匹配 -> 统一效应等位基因方向 -> 处理回文 SNP -> 保存 harmonised 数据
#
# 输出: 02.analysis/opengwas/ldl_chd_harmonised.csv

options(width = 150)
dir.create("02.analysis/opengwas", showWarnings = FALSE)
suppressMessages(library(data.table))

cat("========================================\n")
cat("公开 GWAS 数据 harmonise: LDL-C -> CHD\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 读取暴露工具变量 ----------
cat("[1] 读取 LDL-C 工具变量（68 个，LD clumping 后）\n")
ldl <- fread("02.analysis/ldl_instruments_raw.txt", header = TRUE)
colnames(ldl) <- c("rsid", "A1", "A2", "beta.x", "se.x", "N.x", "P.x", "EAF.x")
cat("    工具变量数:", nrow(ldl), "\n")

# ---------- 2. 读取结局数据 ----------
cat("[2] 读取 CHD 结局数据（CARDIoGRAMplusC4D 2015 完整版）\n")
# 完整文件 792MB 含全基因组; 先按工具变量 rsid 提取再读取, 避免全量载入
ldl_snps <- ldl$rsid
chd_raw <- fread(cmd = paste("zcat -f 00.data/gwas/CHD_CAD2015.txt | awk -F'\\t' 'NR==1 || $1==\"", 
                             paste(ldl_snps, collapse = "\" || $1==\""), "\"'", sep = ""),
                 header = TRUE)
cat("    全基因组行数（完整文件）: ~7.4M; 提取工具变量相关行:", nrow(chd_raw), "\n")
chd <- chd_raw[, .(rsid = markername, A1.c = effect_allele, A2.c = noneffect_allele,
                   beta.y = beta, se.y = se_dgc, P.y = p_dgc,
                   EAF.y = effect_allele_freq)]
chd <- chd[!is.na(beta.y) & !is.na(se.y)]

# ---------- 3. 按 rsid 匹配 ----------
cat("[3] 按 rsid 匹配暴露与结局\n")
m <- merge(ldl, chd, by = "rsid")
cat("    匹配到结局的 SNP:", nrow(m), "/", nrow(ldl), "\n")

# ---------- 4. harmonise: 统一效应等位基因方向 ----------
cat("[4] harmonise（统一效应等位基因方向）\n")
# 暴露 A1 为效应等位基因; 结局 A1.c 为效应等位基因
# 方向一致: 直接保留; 方向相反: beta.y 取负
m[, aligned := ifelse(A1 == A1.c, "same",
               ifelse(A1 == A2.c, "flip", "mismatch"))]
m[aligned == "flip", `:=`(beta.y = -beta.y, EAF.y = 1 - EAF.y, A1.c = A1)]
cat("    方向一致:", sum(m$aligned == "same"),
    "| 翻转后一致:", sum(m$aligned == "flip"),
    "| 无法对齐:", sum(m$aligned == "mismatch"), "\n")
m <- m[aligned != "mismatch"]

# ---------- 5. 检查回文 SNP ----------
cat("[5] 回文 SNP 检查\n")
pal <- m[A1 %in% c("A", "T") & A2 %in% c("A", "T") |
         A1 %in% c("C", "G") & A2 %in% c("C", "G")]
cat("    回文 SNP 数:", nrow(pal), "（无法判断链，需 EAF 辅助或剔除）\n")
# 教学流程: 保留并注明（真实研究中用 EAF 判断链或剔除）
write.csv(pal$rsid, "02.analysis/opengwas/palindromic_snps.txt", row.names = FALSE)

# ---------- 6. 输出 ----------
dat <- m[, .(rsid, A1, A2, beta.exposure = beta.x, se.exposure = se.x,
             pval.exposure = P.x, eaf.exposure = EAF.x,
             beta.outcome = beta.y, se.outcome = se.y, pval.outcome = P.y,
             eaf.outcome = EAF.y, aligned)]
dat$F_stat <- (dat$beta.exposure / dat$se.exposure)^2
cat("\n[6] harmonise 完成\n")
cat("    最终 SNP 数:", nrow(dat), "\n")
cat("    F 统计量: min =", round(min(dat$F_stat), 1),
    ", mean =", round(mean(dat$F_stat), 1), "\n")
write.csv(dat, "02.analysis/opengwas/ldl_chd_harmonised.csv", row.names = FALSE)
cat("    已保存 02.analysis/opengwas/ldl_chd_harmonised.csv\n")
print(dat[, .(rsid, beta.exposure, beta.outcome, F_stat)][1:8], row.names = FALSE)
cat("\nharmonise 完成 ✔\n")
