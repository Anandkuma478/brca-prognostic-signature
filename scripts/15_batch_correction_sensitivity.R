# 15_batch_correction_sensitivity.R
# Sensitivity analysis: PCA before and after removing technical batch, and the frozen models on corrected data.
# Batch is removed with limma::removeBatchEffect while preserving PAM50 subtype. The primary data, models and results are not altered.
# Run from the project root, after 09_metabric_validation.R:
#   Rscript scripts/15_batch_correction_sensitivity.R 2>&1 | tee batch_correction.log

suppressPackageStartupMessages({ library(limma); library(survival); library(data.table); library(ggplot2) })

SERIES <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6")   # fixed categorical order
OTHER  <- "#c3c2b7"

top_genes <- function(m, n = 500) { v <- apply(m, 1, var); names(sort(v, decreasing = TRUE))[seq_len(n)] }
pcs <- function(m, g) {
  p <- prcomp(t(m[g, ]), center = TRUE, scale. = FALSE)
  list(x = p$x[, 1:5], var = round(100 * p$sdev[1:5]^2 / sum(p$sdev^2), 1))
}
adj_r2 <- function(pc, f) {
  ok <- !is.na(f); f <- factor(f[ok])
  if (nlevels(f) < 2) return(NA_real_)
  round(summary(lm(pc[ok] ~ f))$adj.r.squared, 3)
}
lab <- function(x) { x <- as.character(x); x[is.na(x)] <- "Unknown"; factor(x) }
correct <- function(m, bio, batch, batch2 = NULL) removeBatchEffect(m, batch = batch, batch2 = batch2, design = model.matrix(~ lab(bio)))
zscore <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
score <- function(z, s) { k <- s$symbol %in% rownames(z); as.numeric(t(z[s$symbol[k], , drop = FALSE]) %*% s$coef[k]) }

r2_table <- function(cohort, before, after, factors) {
  do.call(rbind, lapply(names(factors), function(n) do.call(rbind, lapply(c("Before", "After"), function(st) {
    p <- if (st == "Before") before else after
    r <- sapply(1:5, function(i) adj_r2(p$x[, i], factors[[n]]))
    data.frame(cohort = cohort, factor = n, stage = st, t(setNames(r, paste0("PC", 1:5, "_adjR2"))), max_adjR2 = max(r))
  }))))
}
evaluate <- function(cohort, model, data, m, sig, time, event, covars) {
  lp <- as.numeric(scale(score(zscore(m), sig)))
  y  <- Surv(time, event)
  co <- concordance(y ~ lp, reverse = TRUE); se <- sqrt(co$var)
  u  <- summary(coxph(y ~ lp))$conf.int[1, ]
  dd <- na.omit(data.frame(t = time, e = event, lp = lp, covars))
  a  <- summary(coxph(Surv(t, e) ~ ., data = dd))
  data.frame(cohort = cohort, model = model, data = data, c_index = co$concordance, c_low = co$concordance - 1.96 * se,
             c_high = co$concordance + 1.96 * se, HR_per_SD = u[1], n_adjusted = nrow(dd), adj_HR = a$conf.int["lp", 1],
             adj_low = a$conf.int["lp", 3], adj_high = a$conf.int["lp", 4], adj_p = a$coefficients["lp", "Pr(>|z|)"], row.names = NULL)
}
panel <- function(before, after, group, title, legend_title, file, h = 3.1) {
  d <- rbind(data.frame(PC1 = before$x[, 1], PC2 = before$x[, 2], Group = group, stage = sprintf("Before correction (PC1 %.1f%%, PC2 %.1f%%)", before$var[1], before$var[2])),
             data.frame(PC1 = after$x[, 1],  PC2 = after$x[, 2],  Group = group, stage = sprintf("After correction (PC1 %.1f%%, PC2 %.1f%%)", after$var[1], after$var[2])))
  d$stage <- factor(d$stage, levels = unique(d$stage))
  d <- d[order(d$Group == "Other plates", decreasing = TRUE), ]
  lv <- setdiff(levels(factor(group)), "Other plates")
  cols <- setNames(c(SERIES[seq_along(lv)], if ("Other plates" %in% group) OTHER), c(lv, if ("Other plates" %in% group) "Other plates"))
  g <- ggplot(d, aes(PC1, PC2, colour = Group)) + geom_point(size = 0.55, alpha = 0.65, stroke = 0) +
    facet_wrap(~ stage, scales = "free") + scale_colour_manual(values = cols, name = legend_title) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) + labs(title = title) +
    theme_bw(base_size = 10) + theme(panel.grid = element_line(colour = "#ecebe7", linewidth = 0.3), strip.background = element_blank(),
                                     strip.text = element_text(colour = "#52514e", size = 9), plot.title = element_text(face = "bold", size = 11),
                                     panel.border = element_rect(colour = "#c3c2b7"))
  ggsave(file, g, width = 7.2, height = h, dpi = 300)
}

models <- readRDS("data/signature_scanb_two_models.rds"); sigs <- list(A = models$sigA, B = models$sigB)
dir_fig <- "figures"

## SCAN-B: sequencing platform
sc <- readRDS("data/scanb_clin.rds"); sx <- readRDS("data/scanb_expr.rds")
stopifnot(identical(rownames(sc), colnames(sx)))
sd <- sx[rowMeans(sx > -3.3) >= 0.5, ]
platform <- sc[["instrument model"]]; sp <- as.character(sc$pam50)
sa <- correct(sd, sp, platform); g <- top_genes(sd)
s0 <- pcs(sd, g); s1 <- pcs(sa, g)
out <- r2_table("SCAN-B", s0, s1, list(Platform = platform, PAM50 = sp))
panel(s0, s1, platform, "SCAN-B", "Platform", file.path(dir_fig, "batch_panel_scanb.png"))
rm(sa, sd, sx); invisible(gc())

## TCGA-BRCA: plate and tissue source site
tc <- readRDS("data/clin_tumor_survival.rds"); tx <- readRDS("data/expr_tumor_vst.rds")
stopifnot(identical(rownames(tc), colnames(tx)))
bc <- rownames(tc); site <- substr(bc, 6, 7); plate <- substr(bc, 22, 25); tp <- as.character(tc$pam50)
ta <- correct(tx, tp, plate, site); g <- top_genes(tx)
t0 <- pcs(tx, g); t1 <- pcs(ta, g)
out <- rbind(out, r2_table("TCGA-BRCA", t0, t1, list(Tissue_source_site = site, Plate = plate, PAM50 = tp)))
top <- names(sort(table(plate), decreasing = TRUE))[1:7]
panel(t0, t1, ifelse(plate %in% top, paste("Plate", plate), "Other plates"), "TCGA-BRCA", "Plate (7 largest)", file.path(dir_fig, "batch_panel_tcga.png"))

gi <- as.data.frame(readRDS("data/gene_info.rds")); sym2id <- setNames(gi$gene_id, gi$gene_name)
tg <- unique(c(sigs$A$symbol, sigs$B$symbol)); id <- sym2id[tg]; keep <- !is.na(id) & id %in% rownames(tx)
sub_t <- function(m) { x <- m[id[keep], ]; rownames(x) <- tg[keep]; x }
cov_t <- data.frame(age = tc$age, stage = tc$stage, pam50 = tc$pam50)
sens <- do.call(rbind, lapply(names(sigs), function(m) rbind(
  evaluate("TCGA-BRCA", m, "Original", sub_t(tx), sigs[[m]], tc$OS_months, tc$OS_event, cov_t),
  evaluate("TCGA-BRCA", m, "Batch-corrected", sub_t(ta), sigs[[m]], tc$OS_months, tc$OS_event, cov_t))))
rm(ta, tx); invisible(gc())

## METABRIC: processing cohort
ex <- fread("data/metabric/raw/data_mrna_illumina_microarray.txt", data.table = FALSE)
L  <- readLines("data/metabric/raw/data_clinical_patient.txt", warn = FALSE)
pa <- read.delim(text = L[!startsWith(L, "#")], quote = "", na.strings = c("", "NA"), check.names = FALSE)
all_m <- as.matrix(ex[, -(1:2)]); rownames(all_m) <- make.unique(ex$Hugo_Symbol)
mm <- all_m[complete.cases(all_m), ]
pp <- pa[match(colnames(mm), pa$PATIENT_ID), c("COHORT", "CLAUDIN_SUBTYPE")]
batch_m <- paste("Cohort", pp$COHORT)
ma <- correct(mm, pp$CLAUDIN_SUBTYPE, batch_m); g <- top_genes(mm)
m0 <- pcs(mm, g); m1 <- pcs(ma, g)
out <- rbind(out, r2_table("METABRIC", m0, m1, list(Processing_cohort = batch_m, PAM50 = pp$CLAUDIN_SUBTYPE)))
panel(m0, m1, batch_m, "METABRIC", "Processing cohort", file.path(dir_fig, "batch_panel_metabric.png"))
rm(ma, mm); invisible(gc())

mb <- readRDS("data/metabric_prepared.rds"); d <- mb$clin
needed <- unique(c(sigs$A$symbol, sigs$B$symbol))
sg <- ex[ex$Hugo_Symbol %in% needed, ]
mat <- as.matrix(sg[, -(1:2)]); rownames(mat) <- sg$Hugo_Symbol
if (anyDuplicated(rownames(mat))) { mat <- mat[order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]; mat <- mat[!duplicated(rownames(mat)), , drop = FALSE] }
mat <- mat[, d$PATIENT_ID, drop = FALSE]
ca  <- correct(mat, as.character(d$pam50), paste("Cohort", pa$COHORT[match(d$PATIENT_ID, pa$PATIENT_ID)]))
cov_m <- d[, c("age", "node", "grade", "er", "pam50", "chemo", "endo")]
sens <- rbind(sens, do.call(rbind, lapply(names(sigs), function(m) rbind(
  evaluate("METABRIC", m, "Original", mat, sigs[[m]], d$time, d$os, cov_m),
  evaluate("METABRIC", m, "Batch-corrected", ca, sigs[[m]], d$time, d$os, cov_m)))))

write.csv(out, "results/batch_before_after.csv", row.names = FALSE)
write.csv(sens, "results/batch_correction_sensitivity.csv", row.names = FALSE)
print(out, row.names = FALSE); cat("\n"); print(format(sens, digits = 3), row.names = FALSE)
cat("\nSaved: results/batch_before_after.csv, results/batch_correction_sensitivity.csv, figures/batch_panel_*.png\n")
