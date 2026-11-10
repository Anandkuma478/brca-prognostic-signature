# 18_mutation_association.R
# Phase 4: somatic mutation frequency by frozen risk group in TCGA-BRCA and METABRIC. No model is refitted.
# Genes tested: mutated (non-synonymous) in at least 5% of sequenced patients in that cohort.
# Tests: Fisher exact (high vs low risk) and logistic regression per SD of score, with and without PAM50 adjustment.
# METABRIC used a targeted panel, so tumour mutation burden is analysed in TCGA-BRCA only.
# Run from the project root with the brca_genomics environment, after 17_tcga_maf_download.R:
#   Rscript scripts/18_mutation_association.R 2>&1 | tee mutation.log

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
NONSYN <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
            "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")
MIN_FREQ <- 0.05

models <- readRDS("data/signature_scanb_two_models.rds"); sigs <- list(A = models$sigA, B = models$sigB)
zscore <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
score  <- function(z, s) { k <- s$symbol %in% rownames(z); as.numeric(t(z[s$symbol[k], , drop = FALSE]) %*% s$coef[k]) }

## Frozen SCAN-B cut-offs
sx <- readRDS("data/scanb_expr.rds")
zs <- zscore(sx[intersect(unique(c(sigs$A$symbol, sigs$B$symbol)), rownames(sx)), ]); rm(sx); invisible(gc())
cuts <- lapply(sigs, function(s) median(score(zs, s)))

analyse <- function(cohort, patients, muts, lp, pam50) {
  mm <- unique(muts[muts$patient %in% patients, c("patient", "gene")])
  fr <- sort(table(mm$gene) / length(patients), decreasing = TRUE)
  genes <- names(fr)[fr >= MIN_FREQ]
  cat(cohort, ": sequenced patients", length(patients), "| genes at >=5%:", length(genes), "\n")
  do.call(rbind, lapply(names(sigs), function(m) {
    high <- lp[[m]] > cuts[[m]]; z <- as.numeric(scale(lp[[m]]))
    r <- do.call(rbind, lapply(genes, function(g) {
      y <- as.numeric(patients %in% mm$patient[mm$gene == g])
      f <- fisher.test(table(factor(high, c(FALSE, TRUE)), factor(y, 0:1)))
      u <- summary(glm(y ~ z, family = binomial))$coefficients["z", ]
      d <- na.omit(data.frame(y = y, z = z, pam50 = pam50))
      a <- summary(glm(y ~ z + pam50, family = binomial, data = d))$coefficients["z", ]
      data.frame(cohort = cohort, model = m, gene = g, n = length(y), n_high = sum(high), mutated_low_pct = 100 * mean(y[!high]),
                 mutated_high_pct = 100 * mean(y[high]), fisher_OR = unname(f$estimate), fisher_p = f$p.value,
                 OR_per_SD = exp(u[1]), OR_per_SD_p = u[4], OR_per_SD_PAM50 = exp(a[1]), OR_per_SD_PAM50_p = a[4], row.names = NULL)
    }))
    r$FDR <- p.adjust(r$fisher_p, "BH"); r$FDR_PAM50 <- p.adjust(r$OR_per_SD_PAM50_p, "BH"); r
  }))
}

## TCGA-BRCA
tc  <- readRDS("data/clin_tumor_survival.rds"); tx <- readRDS("data/expr_tumor_vst.rds")
gi  <- as.data.frame(readRDS("data/gene_info.rds")); sym2id <- setNames(gi$gene_id, gi$gene_name)
tg  <- unique(c(sigs$A$symbol, sigs$B$symbol)); id <- sym2id[tg]; keep <- !is.na(id) & id %in% rownames(tx)
zt  <- zscore(tx[id[keep], ]); rownames(zt) <- tg[keep]; rm(tx); invisible(gc())
maf <- readRDS("data/tcga_maf/tcga_brca_maf.rds")
maf <- maf[substr(maf$Tumor_Sample_Barcode, 14, 15) == "01", ]
maf$patient <- substr(maf$Tumor_Sample_Barcode, 1, 12)
seq_t <- intersect(tc$patient, unique(maf$patient)); i <- match(seq_t, tc$patient)
mt <- maf[maf$Variant_Classification %in% NONSYN, ]; mt$gene <- mt$Hugo_Symbol
lp_t <- lapply(sigs, function(s) score(zt, s)[i])
res <- analyse("TCGA-BRCA", seq_t, mt, lp_t, as.character(tc$pam50)[i])

## Tumour mutation burden (TCGA-BRCA whole exome): non-synonymous variants per patient
n_var <- table(factor(mt$patient, levels = seq_t))
tmb <- do.call(rbind, lapply(names(sigs), function(m) {
  ct <- suppressWarnings(cor.test(as.numeric(n_var), lp_t[[m]], method = "spearman", exact = FALSE))
  d  <- data.frame(y = log10(as.numeric(n_var) + 1), s = as.numeric(scale(lp_t[[m]])), pam50 = as.character(tc$pam50)[i])
  a  <- summary(lm(y ~ s + pam50, data = d))$coefficients["s", ]
  data.frame(model = m, n = length(seq_t), median_variants = median(n_var), rho = unname(ct$estimate), p = ct$p.value,
             beta_log10_PAM50_adjusted = a[1], p_PAM50_adjusted = a[4], row.names = NULL)
}))

## METABRIC (targeted panel)
hdr <- readLines("data/metabric/raw/data_mutations.txt", n = 1)
seqd <- strsplit(sub("^#Sequenced_Samples:\\s*", "", hdr), "\\s+")[[1]]
mu <- fread("data/metabric/raw/data_mutations.txt", skip = 1, data.table = FALSE, quote = "")
mu <- mu[mu$Variant_Classification %in% NONSYN, c("Tumor_Sample_Barcode", "Hugo_Symbol")]
names(mu) <- c("patient", "gene")
mb <- readRDS("data/metabric_prepared.rds"); d <- mb$clin
seq_m <- intersect(d$PATIENT_ID, seqd); j <- match(seq_m, d$PATIENT_ID)
lp_m <- lapply(sigs, function(s) score(mb$z, s)[j])
res <- rbind(res, analyse("METABRIC", seq_m, mu, lp_m, as.character(d$pam50)[j]))

write.csv(res, "results/mutation_association.csv", row.names = FALSE)
write.csv(tmb, "results/tmb_association.csv", row.names = FALSE)

## Figure: genes tested in both cohorts
common <- Reduce(intersect, split(res$gene, res$cohort))
p <- res[res$gene %in% common, ]
ord <- p[p$cohort == "TCGA-BRCA" & p$model == "A", ]; ord <- ord$gene[order((ord$mutated_low_pct + ord$mutated_high_pct) / 2)]
p$gene <- factor(p$gene, levels = ord); p$cohort <- factor(p$cohort, levels = c("TCGA-BRCA", "METABRIC")); p$model <- paste("Model", p$model)
long <- rbind(data.frame(p[, c("cohort", "model", "gene")], pct = p$mutated_low_pct, Group = "Low risk"),
              data.frame(p[, c("cohort", "model", "gene")], pct = p$mutated_high_pct, Group = "High risk"))
long$Group <- factor(long$Group, levels = c("Low risk", "High risk"))
p$label <- ifelse(p$FDR < 0.05, "FDR < 0.05", ""); p$xpos <- pmax(p$mutated_low_pct, p$mutated_high_pct) + 2
g <- ggplot(p, aes(y = gene)) + geom_segment(aes(x = mutated_low_pct, xend = mutated_high_pct, yend = gene), colour = "#c3c2b7", linewidth = 0.7) +
  geom_point(data = long, aes(pct, gene, colour = Group), size = 2.2) +
  geom_text(aes(x = xpos, label = label), hjust = 0, size = 3.1, colour = "#52514e") + facet_grid(model ~ cohort) +
  scale_colour_manual(values = c(`Low risk` = "#2a78d6", `High risk` = "#eb6834"), name = NULL) + scale_x_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  labs(x = "Patients with a non-synonymous mutation (%)", y = NULL) + theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "#ecebe7", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "top", panel.border = element_rect(colour = "#c3c2b7"))
ggsave("figures/mutation_by_risk.png", g, width = 7.5, height = 5 + 0.08 * length(common), dpi = 600, bg = "white")

print(format(res[, c("cohort", "model", "gene", "n", "mutated_low_pct", "mutated_high_pct", "FDR", "OR_per_SD_PAM50", "FDR_PAM50")], digits = 3), row.names = FALSE)
cat("\nTMB (TCGA-BRCA):\n"); print(format(tmb, digits = 3), row.names = FALSE)
cat("\nSaved: results/mutation_association.csv, results/tmb_association.csv, figures/mutation_by_risk.png\n")
