# 20_published_and_random_benchmark.R
# Phase 5: benchmark the frozen models against (i) published signatures, computed with genefu as research-based
# approximations of the clinical assays, and (ii) random signatures of the same size drawn from the candidate genes.
# No model is refitted. Random signatures are weighted by univariate Cox coefficients estimated in SCAN-B.
# Run from the project root with the brca_modules environment, after 09_metabric_validation.R:
#   Rscript scripts/20_published_and_random_benchmark.R 2>&1 | tee benchmark.log

suppressPackageStartupMessages({ library(genefu); library(survival); library(data.table); library(org.Hs.eg.db); library(ggplot2) })
set.seed(2035)
N_BOOT <- 500; N_RANDOM <- 1000

data(sig.oncotypedx); data(sig.gene70); data(sig.ggi); data(pam50)
models <- readRDS("data/signature_scanb_two_models.rds"); sigs <- list(A = models$sigA, B = models$sigB)
cand <- readRDS("data/signature_scanb_model.rds")$candidates$symbol
stopifnot(length(cand) == 1688)

zscore <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
score  <- function(z, s) { k <- s$symbol %in% rownames(z); as.numeric(t(z[s$symbol[k], , drop = FALSE]) %*% s$coef[k]) }
dedupe <- function(m) { m <- m[order(rowMeans(m), decreasing = TRUE), , drop = FALSE]; m[!duplicated(rownames(m)), , drop = FALSE] }
cidx   <- function(y, x) concordance(y ~ x, reverse = TRUE)$concordance

## Published scores. Higher values mean worse prognosis for each (checked against each function's risk call).
published <- function(z) {
  sym <- rownames(z)
  eg  <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, sym, "ENTREZID", "SYMBOL", multiVals = "first"))
  ok  <- !is.na(eg) & !duplicated(eg)
  x   <- t(z[ok, ]); ann <- data.frame(probe = colnames(x), EntrezGene.ID = as.integer(eg[ok]), row.names = colnames(x))
  run <- function(n, f) tryCatch(as.numeric(f()), error = function(e) { cat("  ", n, "failed:", conditionMessage(e), "\n"); rep(NA_real_, nrow(x)) })
  ## oncotypedx returns NA if any of its 21 genes is absent; an absent gene contributes zero, as in the frozen models
  miss <- setdiff(sig.oncotypedx$EntrezGene.ID, ann$EntrezGene.ID)
  if (length(miss)) cat("  Oncotype genes absent (set to zero):", sig.oncotypedx$symbol[match(miss, sig.oncotypedx$EntrezGene.ID)], "\n")
  x_on <- x; ann_on <- ann
  if (length(miss)) {
    id_miss <- paste0("absent_", miss)
    x_on <- cbind(x, matrix(0, nrow(x), length(miss), dimnames = list(NULL, id_miss)))
    ann_on <- rbind(ann, data.frame(probe = id_miss, EntrezGene.ID = as.integer(miss), row.names = id_miss))
  }
  list(`Oncotype DX (21-gene)` = run("oncotype", function() oncotypedx(x_on, ann_on, do.mapping = TRUE, verbose = FALSE)$score),
       `MammaPrint (70-gene)`  = run("gene70", function() gene70(x, ann, do.mapping = TRUE, verbose = FALSE)$score),
       `PAM50 ROR-S`           = run("rorS", function() rorS(x, ann, do.mapping = TRUE, verbose = FALSE)$score),
       `Genomic Grade Index`   = run("ggi", function() ggi(x, ann, do.mapping = TRUE, verbose = FALSE)$score))
}

benchmark <- function(cohort, z, time, event, lp_ours) {
  cat(cohort, "\n")
  pub <- published(z)
  sc  <- c(list(`Model A` = lp_ours$A, `Model B` = lp_ours$B), pub)
  y   <- Surv(time, event); n <- length(time)
  boots <- replicate(N_BOOT, sample(n, replace = TRUE), simplify = FALSE)
  bc <- sapply(sc, function(s) if (anyNA(s)) rep(NA_real_, N_BOOT) else sapply(boots, function(i) cidx(y[i], s[i])))
  out <- do.call(rbind, lapply(names(sc), function(nm) {
    s <- sc[[nm]]
    if (anyNA(s)) return(data.frame(cohort = cohort, signature = nm, c_index = NA, c_low = NA, c_high = NA, diff_A_vs_this = NA, diff_low = NA, diff_high = NA,
                                    spearman_with_A = NA, LR_p_A_added_to_this = NA, LR_p_this_added_to_A = NA))
    d <- bc[, "Model A"] - bc[, nm]
    a <- as.numeric(scale(lp_ours$A)); b <- as.numeric(scale(s))
    lr1 <- if (nm == "Model A") NA else anova(coxph(y ~ b), coxph(y ~ b + a))$`Pr(>|Chi|)`[2]
    lr2 <- if (nm == "Model A") NA else anova(coxph(y ~ a), coxph(y ~ a + b))$`Pr(>|Chi|)`[2]
    data.frame(cohort = cohort, signature = nm, c_index = cidx(y, s), c_low = quantile(bc[, nm], 0.025), c_high = quantile(bc[, nm], 0.975),
               diff_A_vs_this = cidx(y, lp_ours$A) - cidx(y, s), diff_low = quantile(d, 0.025), diff_high = quantile(d, 0.975),
               spearman_with_A = suppressWarnings(cor(lp_ours$A, s, method = "spearman")), LR_p_A_added_to_this = lr1, LR_p_this_added_to_A = lr2, row.names = NULL)
  }))
  out
}

## Cohort data
sc <- readRDS("data/scanb_clin.rds"); sx <- readRDS("data/scanb_expr.rds"); cv <- readRDS("data/nested_cv_predictions.rds")
zs <- zscore(sx[rowMeans(sx > -3.3) >= 0.5, ]); rm(sx); invisible(gc())
tc <- readRDS("data/clin_tumor_survival.rds"); tx <- readRDS("data/expr_tumor_vst.rds")
gi <- as.data.frame(readRDS("data/gene_info.rds")); nm <- gi$gene_name[match(rownames(tx), gi$gene_id)]
tx <- tx[!is.na(nm) & nm != "", ]; rownames(tx) <- nm[!is.na(nm) & nm != ""]; zt <- zscore(dedupe(tx)); rm(tx); invisible(gc())
## CD68 (an Oncotype DX gene) was removed from the filtered TCGA matrix; restore it from unfiltered TPM for this benchmark only
if (!"CD68" %in% rownames(zt)) {
  se <- readRDS("data/tcga_brca_se.rds"); i <- which(SummarizedExperiment::rowData(se)$gene_name == "CD68")
  stopifnot(length(i) == 1, all(colnames(zt) %in% colnames(se)))
  cd <- log2(SummarizedExperiment::assay(se, "tpm_unstrand")[i, colnames(zt)] + 1)
  zt <- rbind(zt, CD68 = as.numeric(scale(cd))); rm(se); invisible(gc())
  cat("CD68 restored for TCGA-BRCA from unfiltered log2(TPM + 1)\n")
}
ex <- fread("data/metabric/raw/data_mrna_illumina_microarray.txt", data.table = FALSE)
mb <- readRDS("data/metabric_prepared.rds"); d <- mb$clin
mx <- as.matrix(ex[, d$PATIENT_ID]); rownames(mx) <- ex$Hugo_Symbol; rm(ex); invisible(gc())
mx <- mx[!is.na(rownames(mx)) & rownames(mx) != "", ]
mx <- mx[order(rowMeans(mx, na.rm = TRUE), decreasing = TRUE), ]; mx <- mx[!duplicated(rownames(mx)), ]
zm <- zscore(mx); rm(mx); invisible(gc())

res <- rbind(
  benchmark("SCAN-B (out-of-fold)", zs, sc$OS_months, sc$OS_event, list(A = rowMeans(cv$lpA), B = rowMeans(cv$lpB))),
  benchmark("TCGA-BRCA", zt, tc$OS_months, tc$OS_event, lapply(sigs, function(s) score(zt, s))),
  benchmark("METABRIC", zm, d$time, d$os, lapply(sigs, function(s) score(zm, s))))
write.csv(res, "results/published_signature_comparison.csv", row.names = FALSE)

## Gene overlap with published signatures
lists <- list(`Oncotype DX (21-gene)` = sig.oncotypedx$symbol[sig.oncotypedx$group != "reference"],
              `MammaPrint (70-gene)` = sig.gene70$HUGO.gene.symbol, `PAM50` = rownames(pam50$centroids),
              `Genomic Grade Index` = sig.ggi$HUGO.gene.symbol)
ov <- do.call(rbind, lapply(names(lists), function(n) do.call(rbind, lapply(names(sigs), function(m) {
  l <- unique(na.omit(lists[[n]])); o <- intersect(sigs[[m]]$symbol, l)
  data.frame(published = n, genes_in_published = length(l), model = m, overlap = length(o), genes = paste(o, collapse = ", "))
}))))
write.csv(ov, "results/signature_gene_overlap.csv", row.names = FALSE)

## Random signatures: candidate genes measured in every cohort, weighted by SCAN-B univariate Cox coefficients
pool <- Reduce(intersect, list(cand, rownames(zs), rownames(zt), rownames(zm)))
cat("Random-signature pool:", length(pool), "genes\n")
ys <- Surv(sc$OS_months, sc$OS_event)
ub <- sapply(pool, function(g) coef(coxph(ys ~ zs[g, ])))
names(ub) <- pool
rnd <- do.call(rbind, lapply(names(sigs), function(m) {
  k <- nrow(sigs[[m]])
  sets <- replicate(N_RANDOM, sample(pool, k), simplify = FALSE)
  do.call(rbind, lapply(list(list("TCGA-BRCA", zt, tc$OS_months, tc$OS_event), list("METABRIC", zm, d$time, d$os)), function(co) {
    y <- Surv(co[[3]], co[[4]]); z <- co[[2]]
    obs <- cidx(y, score(z, sigs[[m]]))
    rc  <- sapply(sets, function(g) cidx(y, as.numeric(crossprod(z[g, , drop = FALSE], ub[g]))))
    data.frame(cohort = co[[1]], model = m, genes = k, n_random = N_RANDOM, observed_c = obs, random_median = median(rc), random_p2.5 = unname(quantile(rc, 0.025)),
               random_p97.5 = unname(quantile(rc, 0.975)), share_random_at_least_observed = mean(rc >= obs), empirical_p = (sum(rc >= obs) + 1) / (N_RANDOM + 1))
  }))
}))
write.csv(rnd, "results/random_signature_benchmark.csv", row.names = FALSE)

## Figure
f <- res[!is.na(res$c_index), ]
f$cohort <- factor(f$cohort, levels = c("SCAN-B (out-of-fold)", "TCGA-BRCA", "METABRIC"))
f$signature <- factor(f$signature, levels = rev(c("Model A", "Model B", "Oncotype DX (21-gene)", "MammaPrint (70-gene)", "PAM50 ROR-S", "Genomic Grade Index")))
f$Group <- ifelse(grepl("^Model", f$signature), "This study", "Published (research-based)")
g <- ggplot(f, aes(c_index, signature, colour = Group)) + geom_vline(xintercept = 0.5, linetype = 2, colour = "#c3c2b7") +
  geom_segment(aes(x = c_low, xend = c_high, yend = signature), linewidth = 0.7) + geom_point(size = 2.2) +
  geom_text(aes(x = 0.81, label = sprintf("%.3f", c_index)), colour = "#52514e", size = 2.8, hjust = 1, show.legend = FALSE) +
  facet_wrap(~ cohort) + scale_colour_manual(values = c(`This study` = "#2a78d6", `Published (research-based)` = "#eb6834"), name = NULL) +
  scale_x_continuous(limits = c(0.45, 0.82), breaks = c(0.5, 0.6, 0.7, 0.8)) + labs(x = "C-index for overall survival (95% bootstrap CI)", y = NULL) + theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "#ecebe7", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "top", panel.border = element_rect(colour = "#c3c2b7"))
ggsave("figures/published_signature_forest.png", g, width = 8, height = 3.3, dpi = 600, bg = "white")

options(width = 170)
print(format(res, digits = 3), row.names = FALSE); cat("\n"); print(ov, row.names = FALSE); cat("\n"); print(format(rnd, digits = 3), row.names = FALSE)
cat("\nSaved: results/published_signature_comparison.csv, signature_gene_overlap.csv, random_signature_benchmark.csv, figures/published_signature_forest.png\n")
