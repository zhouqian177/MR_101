#!/usr/bin/env Rscript
# 36_eqtl_mapping.R — eQTL 定位：基因表达 -> 疾病（转录组 MR + coloc 共定位）
#
# 背景: eQTL 定位（Transcriptome-wide MR / colocalization）用基因表达 eQTL
#       作为工具变量, 评估"基因表达水平"对疾病的因果效应, 并可与疾病 GWAS
#       做共定位（coloc）区分"共享因果变异"与"LD 假关联"。
#
# 实例: HMGCR（HMG-CoA 还原酶, 他汀类药物的靶点基因）基因表达 -> CHD
#   - 暴露: eqtl-a-ENSG00000140464 (HMGCR 基因表达 eQTL, OpenGWAS)
#   - 结局: ieu-a-7 (CHD)
#   - cis 工具变量: rs10851868 (P=3.2e-145), rs12904134 (P=9.9e-10)
#
# 方法: 1) 基因表达 -> CHD 的 MR（IVW/单工具 Wald ratio）
#       2) coloc 共定位（简化演示, 使用 eQTL 与结局的共享 SNP）
#
# 输出: 02.analysis/multiomics/eqtl_*.csv + 04.figures/eqtl_*.png

options(width = 150)
dir.create("02.analysis/multiomics", showWarnings = FALSE)
suppressMessages(library(TwoSampleMR))
suppressMessages(library(coloc))

cat("========================================\n")
cat("eQTL 定位: HMGCR 基因表达 -> CHD\n")
cat("时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# ---------- 1. 数据 ----------
cat("[1] 数据（OpenGWAS 在线）\n")
eqtl <- read.table("02.analysis/multiomics/eqtl_HMGCR.csv", header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
chd  <- read.table("02.analysis/multiomics/outcome_HMGCR_CHD.csv", header = TRUE,
                   sep = "\t", stringsAsFactors = FALSE)
cat("    HMGCR eQTL 工具变量:", nrow(eqtl), "个\n")
cat("    CHD 关联:", nrow(chd), "条\n")

# ---------- 2. 基因表达 -> CHD MR ----------
cat("\n[2] 转录组 MR: HMGCR 表达 -> CHD\n")
m <- merge(eqtl, chd, by = "rsid")
# 统一效应等位基因方向（merge 后列名 ea.x/ea.y）
m$flip <- m$ea.x != m$ea.y
m$beta_y <- ifelse(m$flip, -m$beta.y, m$beta.y)
m$se_y <- m$se.y
m$beta_x <- m$beta.x; m$se_x <- m$se.x
# IVW 汇总（2 个工具）
w <- 1 / m$se_y^2
ivw_b <- sum(m$beta_x * m$beta_y * w) / sum(m$beta_x^2 * w)
ivw_se <- sqrt(1 / sum(m$beta_x^2 * w))
ivw_p <- 2 * pnorm(-abs(ivw_b / ivw_se))
cat(sprintf("    IVW: HMGCR 表达每升高 1 SD, CHD logOR = %.4f (SE=%.4f, P=%.3g)\n",
            ivw_b, ivw_se, ivw_p))
cat(sprintf("    OR = %.3f (%.3f-%.3f)\n", exp(ivw_b),
            exp(ivw_b - 1.96 * ivw_se), exp(ivw_b + 1.96 * ivw_se)))
# 单工具 Wald ratio
for (k in 1:nrow(m)) {
  b <- m$beta_y[k] / m$beta_x[k]; se <- m$se_y[k] / abs(m$beta_x[k])
  cat(sprintf("    工具 %s: Wald ratio = %.4f (SE=%.4f, P=%.3g)\n",
              m$rsid[k], b, se, 2 * pnorm(-abs(b / se))))
}
res <- data.frame(方法 = c("IVW", paste0("Wald_", m$rsid)),
                  beta = c(ivw_b, m$beta_y / m$beta_x),
                  SE = c(ivw_se, m$se_y / abs(m$beta_x)),
                  P = c(ivw_p, 2 * pnorm(-abs((m$beta_y / m$beta_x) / (m$se_y / abs(m$beta_x))))))
write.csv(res, "02.analysis/multiomics/eqtl_HMGCR_mr.csv", row.names = FALSE)

# ---------- 3. coloc 共定位（简化演示） ----------
cat("\n[3] coloc 共定位（eQTL 与 CHD 共享因果变异检验）\n")
# 用 eQTL 区域工具变量 + CHD 关联做 coloc.abf
# 需要 beta/varbeta/MAF/type/N；MAF 需为严格 (0,1) 数值
maf1 <- as.numeric(m$eaf.x); maf2 <- as.numeric(m$eaf.y)
maf1 <- pmin(pmax(maf1, 1e-4), 1 - 1e-4)
maf2 <- pmin(pmax(maf2, 1e-4), 1 - 1e-4)
d1 <- list(beta = as.numeric(m$beta_x), varbeta = as.numeric(m$se_x)^2,
           MAF = maf1, type = "quant", N = 200, snp = m$rsid)      # eQTL
d2 <- list(beta = as.numeric(m$beta_y), varbeta = as.numeric(m$se_y)^2,
           MAF = maf2, type = "cc", s = 60801 / (60801 + 123504), # CHD 病例比例
           N = 184305, snp = m$rsid)
col <- tryCatch(coloc.abf(d1, d2), error = function(e) {
  cat("    coloc 失败:", conditionMessage(e), "\n"); NULL
})
if (!is.null(col)) {
  pp <- col$summary
  cat(sprintf("    PP.H0=%.3f PP.H1=%.3f PP.H2=%.3f PP.H3=%.3f PP.H4=%.3f\n",
              pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"],
              pp["PP.H3.abf"], pp["PP.H4.abf"]))
  cat("    PP.H4 > 0.8 支持 eQTL 与疾病共享因果变异（真共定位）\n")
}

# ---------- 4. 解读 ----------
cat("\n[4] 解读\n")
cat("    - HMGCR 是他汀靶点基因; 他汀降低 HMGCR 活性 -> 降低 LDL-C 与 CHD 风险\n")
cat("    - 若 eQTL 定位显示 HMGCR 表达增加 CHD 风险, 支持该靶点的药物重定位价值\n")
cat("\neQTL 定位完成 ✔\n")
