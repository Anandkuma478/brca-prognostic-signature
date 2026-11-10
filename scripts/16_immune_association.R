# 16_immune_association.R
# Phase 4: association of the frozen risk scores with immune features in SCAN-B, TCGA-BRCA and METABRIC.
# Features: ssGSEA of immune-related Hallmark sets, a cytolytic score (GZMA, PRF1), CD8A and immune checkpoint genes.
# Run from the project root with the brca_modules environment, after 09_metabric_validation.R:
#   Rscript scripts/16_immune_association.R 2>&1 | tee immune.log

suppressPackageStartupMessages({ library(GSVA); library(data.table); library(msigdbr); library(ggplot2); library(BiocParallel) })
set.seed(2034)
PAR <- MulticoreParam(8)

SETS <- paste0("HALLMARK_", c("INTERFERON_GAMMA_RESPONSE", "INTERFERON_ALPHA_RESPONSE", "INFLAMMATORY_RESPONSE", "ALLOGRAFT_REJECTION",
                              "IL2_STAT5_SIGNALING", "IL6_JAK_STAT3_SIGNALING", "TNFA_SIGNALING_VIA_NFKB", "COMPLEMENT"))
CHECK <- c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1")

h  <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "H"))
h  <- h[h$gs_name %in% SETS, c("gs_name", "gene_symbol")]
gs <- lapply(split(h$gene_symbol, h$gs_name), unique)
models <- readRDS("data/signature_scanb_two_models.rds"); sigs <- list(A = models$sigA, B = models$sigB)
cat("Model genes also in the immune Hallmark sets: A", sum(sigs$A$symbol %in% unlist(gs)), "of", nrow(sigs$A),
    "| B", sum(sigs$B$symbol %in% unlist(gs)), "of", nrow(sigs$B), "\n")

zscore <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
score  <- function(z, s) { k <- s$symbol %in% rownames(z); as.numeric(t(z[s$symbol[k], , drop = FALSE]) %*% s$coef[k]) }
dedupe <- function(m) { m <- m[order(rowMeans(m), decreasing = TRUE), , drop = FALSE]; m[!duplicated(rownames(m)), , drop = FALSE] }

features <- function(m) {
  ss <- gsva(ssgseaParam(m, gs), BPPARAM = PAR, verbose = FALSE)
  rownames(ss) <- sub("HALLMARK_", "", rownames(ss))
  z  <- zscore(m[intersect(c("GZMA", "PRF1", "CD8A", CHECK), rownames(m)), , drop = FALSE])
  stopifnot(all(c("GZMA", "PRF1") %in% rownames(z)))
  rbind(ss, Cytolytic_score = colMeans(z[c("GZMA", "PRF1"), ]), z[setdiff(rownames(z), c("GZMA", "PRF1")), , drop = FALSE])
}
analyse <- function(cohort, f, lp, pam50) {
  out <- do.call(rbind, lapply(names(lp), function(m) {
    r <- do.call(rbind, lapply(rownames(f), function(x) {
      ct <- suppressWarnings(cor.test(f[x, ], lp[[m]], method = "spearman", exact = FALSE))
      d  <- data.frame(y = as.numeric(scale(f[x, ])), s = as.numeric(scale(lp[[m]])), pam50 = pam50)
      a  <- summary(lm(y ~ s + pam50, data = d))$coefficients["s", ]
      data.frame(cohort = cohort, model = m, feature = x, rho = unname(ct$estimate), p = ct$p.value,
                 beta_PAM50_adjusted = a[1], p_PAM50_adjusted = a[4], row.names = NULL)
    }))
    r$FDR <- p.adjust(r$p, "BH"); r$FDR_PAM50_adjusted <- p.adjust(r$p_PAM50_adjusted, "BH"); r
  }))
  cat(cohort, "done\n"); out
}

## SCAN-B: out-of-fold scores
sc <- readRDS("data/scanb_clin.rds"); sx <- readRDS("data/scanb_expr.rds"); cv <- readRDS("data/nested_cv_predictions.rds")
stopifnot(identical(rownames(sc), colnames(sx)))
res <- analyse("SCAN-B", features(sx[rowMeans(sx > -3.3) >= 0.5, ]), list(A = rowMeans(cv$lpA), B = rowMeans(cv$lpB)), as.character(sc$pam50))
rm(sx); invisible(gc())

## TCGA-BRCA: frozen models
tc <- readRDS("data/clin_tumor_survival.rds"); tx <- readRDS("data/expr_tumor_vst.rds")
gi <- as.data.frame(readRDS("data/gene_info.rds")); nm <- gi$gene_name[match(rownames(tx), gi$gene_id)]
tx <- tx[!is.na(nm) & nm != "", ]; rownames(tx) <- nm[!is.na(nm) & nm != ""]; tx <- dedupe(tx)
res <- rbind(res, analyse("TCGA-BRCA", features(tx), lapply(sigs, function(s) score(zscore(tx[intersect(s$symbol, rownames(tx)), ]), s)), as.character(tc$pam50)))
rm(tx); invisible(gc())

## METABRIC: frozen models
ex <- fread("data/metabric/raw/data_mrna_illumina_microarray.txt", data.table = FALSE)
mb <- readRDS("data/metabric_prepared.rds"); d <- mb$clin
mx <- as.matrix(ex[, d$PATIENT_ID]); rownames(mx) <- ex$Hugo_Symbol
mx <- mx[!is.na(rownames(mx)) & rownames(mx) != "" & complete.cases(mx), ]; mx <- dedupe(mx)
res <- rbind(res, analyse("METABRIC", features(mx), lapply(sigs, function(s) score(mb$z, s)), as.character(d$pam50)))

write.csv(res, "results/immune_association.csv", row.names = FALSE)

## Figure
pretty <- function(x) { x <- gsub("_", " ", x); x <- sub("Cytolytic score", "Cytolytic score (GZMA, PRF1)", x); x }
order_f <- c("Cytolytic_score", "CD8A", CHECK, sub("HALLMARK_", "", SETS))
res$feature_lab <- factor(pretty(res$feature), levels = rev(pretty(order_f)))
res$cohort <- factor(res$cohort, levels = c("SCAN-B", "TCGA-BRCA", "METABRIC"))
res$panel <- paste("Model", res$model)
g <- ggplot(res, aes(cohort, feature_lab, fill = rho)) + geom_tile(colour = "#fcfcfb", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.4, colour = "#0b0b0b") + facet_wrap(~ panel) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-0.8, 0.8), oob = scales::squish, name = "Spearman rho") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"), axis.text.y = element_text(colour = "#0b0b0b"))
ggsave("figures/immune_association_heatmap.png", g, width = 8, height = 5.6, dpi = 600, bg = "white")

sm <- aggregate(cbind(rho, beta_PAM50_adjusted) ~ cohort + model, res, median)
print(format(res[, c("cohort", "model", "feature", "rho", "FDR", "beta_PAM50_adjusted", "FDR_PAM50_adjusted")], digits = 3), row.names = FALSE)
cat("\nMedian over features:\n"); print(format(sm, digits = 3), row.names = FALSE)
cat("\nSaved: results/immune_association.csv, figures/immune_association_heatmap.png\n")
