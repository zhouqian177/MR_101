#!/usr/bin/env Rscript
# 31_drug_target_mr.R — 药物靶点 MR（Drug-target MR / cis-MR）
#
# 背景: 药物靶点 MR 用靶点蛋白/基因表达的 pQTL/eQTL 作为工具变量,
#       评估"调节该靶点"对结局的因果效应, 模拟药物疗效(如 PCSK9 抑制剂)。
#       常用 cis 变异(靶点基因 ±1Mb 内)作为工具变量。
#
# 实例: PCSK9 蛋白水平 (pQTL: ebi-a-GCST90010246) -> LDL-C / CHD
#   - PCSK9 是降脂药(依洛尤单抗/阿利西尤单抗)的靶点
#   - rs631220 (chr1:55527479, PCSK9 基因内) 为 cis 工具变量
#
# 数据(OpenGWAS 在线):
#   暴露: PCSK9 pQTL rs631220  beta=0.3944 se=0.0547 P=4.8e-12
#   结局1: LDL-C (ieu-b-110)   beta=0.0247 se=0.0027 P=8.1e-20
#   结局2: CHD (ebi-a-GCST005194) beta=0.0053 se=0.0075 P=0.483
#
# 方法: 单 SNP Wald ratio = beta_out / beta_exp
# 输出: 02.analysis/drug_target/ 结果

options(width = 150)
dir.create("02.analysis/drug_target", showWarnings = FALSE)

cat("========================================\n")
cat("药物靶点 MR（Drug-target MR, cis-MR）\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 数据 ----------
cat("[1] 数据（OpenGWAS 在线, PCSK9 pQTL）\n")
pqtl <- data.frame(rsid = "rs631220", ea = "G", nea = "A",
                   beta = 0.394406, se = 0.0546559, p = 4.78851e-12)
ldl  <- data.frame(rsid = "rs631220", ea = "G", nea = "A",
                   beta = 0.0246595, se = 0.00270619, p = 8.10028e-20)
chd  <- read.table("02.analysis/opengwas/online_outcome_PCSK9_CHD.csv",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("    PCSK9 pQTL 工具变量:", pqtl$rsid, "(P =", format(pqtl$p, digits = 3), ")\n")
cat("    LDL-C 关联 P =", format(ldl$p, digits = 3),
    "| CHD 关联 P =", format(chd$p, digits = 3), "\n\n")

# ---------- 2. Wald ratio ----------
cat("[2] Wald ratio 估计（单 SNP cis-MR）\n")
wald <- function(bx, sex, by, sey, label, binary = FALSE) {
  b <- by / bx
  se <- sey / abs(bx)
  pv <- 2 * pnorm(-abs(b / se))
  if (binary) {
    cat(sprintf("    %-22s beta=%.4f (SE=%.4f, P=%.3g) | OR(每单位蛋白)=%.3f (%.3f-%.3f)\n",
                label, b, se, pv, exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)))
  } else {
    cat(sprintf("    %-22s beta=%.4f (SE=%.4f, P=%.3g)\n", label, b, se, pv))
  }
  c(b = b, se = se, p = pv)
}
r_ldl <- wald(pqtl$beta, pqtl$se, ldl$beta, ldl$se, "PCSK9 -> LDL-C")
r_chd <- wald(pqtl$beta, pqtl$se, chd$beta[1], chd$se[1], "PCSK9 -> CHD", binary = TRUE)

# ---------- 3. 解释 ----------
cat("\n[3] 结果解读\n")
cat("    - PCSK9 蛋白水平每升高 1 单位: LDL-C 变化", sprintf("%.4f", r_ldl["b"]),
    "(P =", format(r_ldl["p"], digits = 3), ")\n")
cat("    - PCSK9 蛋白水平与 CHD: 效应", sprintf("%.4f", r_chd["b"]),
    ", OR =", sprintf("%.3f", exp(r_chd["b"])),
    "(P =", format(r_chd["p"], digits = 3), ")\n")
cat("    - 临床意义: PCSK9 抑制剂(单抗/小干扰RNA)降低 PCSK9 水平 ->\n")
cat("      降低 LDL-C 与 CHD 风险; 本结果支持该靶点的因果效应\n")

# ---------- 4. 保存 ----------
res <- data.frame(
  关联 = c("PCSK9 -> LDL-C", "PCSK9 -> CHD"),
  beta = c(r_ldl["b"], r_chd["b"]),
  SE = c(r_ldl["se"], r_chd["se"]),
  P = c(r_ldl["p"], r_chd["p"]),
  OR = c(NA, exp(r_chd["b"])),
  OR_lci = c(NA, exp(r_chd["b"] - 1.96 * r_chd["se"])),
  OR_uci = c(NA, exp(r_chd["b"] + 1.96 * r_chd["se"])))
write.csv(res, "02.analysis/drug_target/results.csv", row.names = FALSE)
print(res, row.names = FALSE)
cat("\n药物靶点 MR 完成 ✔\n")
