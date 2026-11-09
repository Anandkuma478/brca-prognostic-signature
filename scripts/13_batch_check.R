# 13_batch_check.R
# Quality check: how much expression variation tracks technical batch, compared with biology (PAM50).
# SCAN-B: sequencing platform. TCGA-BRCA tumours: tissue source site and plate from the sample barcode.
# Reports adjusted R-squared of the first five principal components (500 most variable genes). Data are not altered.
# Run from the project root:  Rscript scripts/13_batch_check.R 2>&1 | tee batch_check.log

suppressPackageStartupMessages(library(ggplot2))

pcs <- function(m, n_top = 500) {
  v <- apply(m, 1, var); top <- names(sort(v, decreasing = TRUE))[seq_len(n_top)]
  p <- prcomp(t(m[top, ]), center = TRUE, scale. = FALSE)
  list(x = p$x[, 1:5], var = round(100 * p$sdev[1:5]^2 / sum(p$sdev^2), 1))
}
adj_r2 <- function(pc, f) {
  ok <- !is.na(f); f <- factor(f[ok])
  if (nlevels(f) < 2) return(NA_real_)
  round(summary(lm(pc[ok] ~ f))$adj.r.squared, 3)
}
report <- function(cohort, p, factors) {
  do.call(rbind, lapply(names(factors), function(n) {
    data.frame(cohort = cohort, factor = n, levels = length(unique(na.omit(factors[[n]]))),
               t(setNames(sapply(1:5, function(i) adj_r2(p$x[, i], factors[[n]])), paste0("PC", 1:5, "_adjR2"))))
  }))
}

sc <- readRDS("data/scanb_clin.rds"); sx <- readRDS("data/scanb_expr.rds")
stopifnot(identical(rownames(sc), colnames(sx)))
det <- rowMeans(sx > -3.3) >= 0.5   # without this, the most variable genes are mostly at the expression floor
cat("SCAN-B genes detected in at least 50% of samples:", sum(det), "\n")
p_s <- pcs(sx[det, ])
platform <- sc[["instrument model"]]
cat("SCAN-B platform counts:\n"); print(table(platform))
cat("SCAN-B variance explained, PC1-5 (%):", p_s$var, "\n")
out_s <- report("SCAN-B", p_s, list(Platform = platform, PAM50 = as.character(sc$pam50)))

tc <- readRDS("data/clin_tumor_survival.rds"); tx <- readRDS("data/expr_tumor_vst.rds")
stopifnot(identical(rownames(tc), colnames(tx)))
bc <- rownames(tc)
p_t <- pcs(tx)
cat("TCGA variance explained, PC1-5 (%):", p_t$var, "\n")
out_t <- report("TCGA-BRCA", p_t, list(Tissue_source_site = substr(bc, 6, 7), Plate = substr(bc, 22, 25),
                                       PAM50 = as.character(tc$pam50)))
## METABRIC: processing cohort. Expression and cohort labels only; survival is not read.
ex <- data.table::fread("data/metabric/raw/data_mrna_illumina_microarray.txt", data.table = FALSE)
mm <- as.matrix(ex[, -(1:2)]); rownames(mm) <- make.unique(ex$Hugo_Symbol)
mm <- mm[complete.cases(mm), ]
L  <- readLines("data/metabric/raw/data_clinical_patient.txt", warn = FALSE)
pa <- read.delim(text = L[!startsWith(L, "#")], quote = "", na.strings = c("", "NA"), check.names = FALSE)
pa <- pa[match(colnames(mm), pa$PATIENT_ID), c("COHORT", "CLAUDIN_SUBTYPE")]
p_m <- pcs(mm)
cat("METABRIC variance explained, PC1-5 (%):", p_m$var, "\n")
out_m <- report("METABRIC", p_m, list(Processing_cohort = as.character(pa$COHORT), PAM50 = pa$CLAUDIN_SUBTYPE))
out <- rbind(out_s, out_t, out_m)
write.csv(out, "results/batch_check.csv", row.names = FALSE)
print(out, row.names = FALSE)

d <- data.frame(PC1 = p_s$x[, 1], PC2 = p_s$x[, 2], Platform = platform)
g <- ggplot(d, aes(PC1, PC2, colour = Platform)) + geom_point(size = 0.7, alpha = 0.6) +
  labs(x = sprintf("PC1: %.1f%%", p_s$var[1]), y = sprintf("PC2: %.1f%%", p_s$var[2]), title = "SCAN-B by sequencing platform") +
  theme_bw(base_size = 11)
ggsave("figures/batch_pca_scanb.png", g, width = 6, height = 5, dpi = 300)
cat("\nSaved: results/batch_check.csv, figures/batch_pca_scanb.png\n")
