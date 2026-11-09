# 08_wgcna.R
# Phase 2: WGCNA co-expression modules in TCGA-BRCA tumours, module-trait association,
# and overlap of each module with the tumour-vs-normal DEGs.
# Run from the project root:
#   Rscript scripts/08_wgcna.R 2>&1 | tee wgcna.log
# Optional arguments: number of top-variance genes (default 5000), output directory (default ".").

args   <- commandArgs(trailingOnly = TRUE)
n_top  <- if (length(args) >= 1) as.integer(args[1]) else 5000L
outdir <- if (length(args) >= 2) args[2] else "."
resdir <- file.path(outdir, "results"); figdir <- file.path(outdir, "figures")
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 8)
set.seed(2029)
FALLBACK_POWER <- 12   # WGCNA FAQ value for a signed network with more than 40 samples

expr <- readRDS("data/expr_tumor_vst.rds")
clin <- readRDS("data/clin_tumor_survival.rds")
gi   <- as.data.frame(readRDS("data/gene_info.rds"))
stopifnot(identical(rownames(clin), colnames(expr)))
cat("Tumours:", ncol(expr), " Genes in matrix:", nrow(expr), "\n")

## 1. Protein-coding genes, top variance
pc    <- gi$gene_id[gi$gene_type == "protein_coding"]
genes <- intersect(rownames(expr), pc)
v     <- apply(expr[genes, ], 1, var)
top   <- names(sort(v, decreasing = TRUE))[seq_len(min(n_top, length(v)))]
dat   <- t(expr[top, ])
gsg   <- goodSamplesGenes(dat, verbose = 0)
stopifnot(gsg$allOK)
cat("Protein-coding genes:", length(genes), " Used for WGCNA:", ncol(dat), "\n")

## 2. Sample clustering (reported, no samples are removed automatically)
tree <- hclust(dist(dat), method = "average")
png(file.path(figdir, "wgcna_sample_tree.png"), width = 2400, height = 1200, res = 200)
plot(tree, labels = FALSE, main = "Sample clustering (TCGA-BRCA tumours)", xlab = "", sub = "")
dev.off()

## 3. Soft-threshold power
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(dat, powerVector = powers, networkType = "signed", RsquaredCut = 0.85, verbose = 0)
fit <- data.frame(power = sft$fitIndices$Power,
                  scale_free_R2 = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
                  mean_connectivity = sft$fitIndices$mean.k.)
write.csv(fit, file.path(resdir, "wgcna_soft_threshold.csv"), row.names = FALSE)
power <- if (is.na(sft$powerEstimate)) FALLBACK_POWER else sft$powerEstimate
cat("Power estimate:", sft$powerEstimate, " Power used:", power,
    " Scale-free R2 at power:", round(fit$scale_free_R2[fit$power == power], 3), "\n")

png(file.path(figdir, "wgcna_soft_threshold.png"), width = 2400, height = 1200, res = 200)
par(mfrow = c(1, 2))
plot(fit$power, fit$scale_free_R2, type = "n", xlab = "Soft threshold power", ylab = "Scale-free fit (signed R2)")
text(fit$power, fit$scale_free_R2, labels = fit$power, col = "red"); abline(h = 0.85, col = "red", lty = 2)
plot(fit$power, fit$mean_connectivity, type = "n", xlab = "Soft threshold power", ylab = "Mean connectivity")
text(fit$power, fit$mean_connectivity, labels = fit$power, col = "red")
dev.off()

## 4. Network and modules
cor <- WGCNA::cor
net <- blockwiseModules(dat, power = power, networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, mergeCutHeight = 0.25, maxBlockSize = ncol(dat) + 1,
                        numericLabels = TRUE, saveTOMs = FALSE, randomSeed = 2029, verbose = 0)
cor <- stats::cor
colors <- labels2colors(net$colors)
names(colors) <- colnames(dat)
cat("Modules (excluding grey):", length(setdiff(unique(colors), "grey")),
    " Unassigned (grey) genes:", sum(colors == "grey"), "\n")
print(sort(table(colors), decreasing = TRUE))

png(file.path(figdir, "wgcna_dendrogram.png"), width = 2400, height = 1200, res = 200)
plotDendroAndColors(net$dendrograms[[1]], colors[net$blockGenes[[1]]], "Module", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules")
dev.off()

## 5. Module-trait association
MEs <- orderMEs(moduleEigengenes(dat, colors)$eigengenes)
pam <- as.character(clin$pam50)
traits <- data.frame(Age = clin$age, Stage = as.numeric(clin$stage),
                     Basal = as.numeric(pam == "Basal"), Her2 = as.numeric(pam == "Her2"),
                     LumA = as.numeric(pam == "LumA"), LumB = as.numeric(pam == "LumB"),
                     Normal_like = as.numeric(pam == "Normal"), row.names = rownames(clin))
cp  <- corAndPvalue(MEs, traits, use = "p")
fdr <- matrix(p.adjust(cp$p, "BH"), nrow(cp$p), dimnames = dimnames(cp$p))
long <- do.call(rbind, lapply(colnames(traits), function(t)
  data.frame(module = sub("^ME", "", rownames(cp$cor)), trait = t, cor = cp$cor[, t], p = cp$p[, t], FDR = fdr[, t])))
write.csv(long[order(long$FDR), ], file.path(resdir, "wgcna_module_trait.csv"), row.names = FALSE)

png(file.path(figdir, "wgcna_module_trait.png"), width = 2000, height = 2000, res = 220)
par(mar = c(6, 9, 3, 2))
labeledHeatmap(Matrix = cp$cor, xLabels = colnames(traits), yLabels = colnames(MEs), ySymbols = colnames(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50), zlim = c(-1, 1),
               textMatrix = paste0(signif(cp$cor, 2), "\n(", signif(fdr, 1), ")"),
               setStdMargins = FALSE, cex.text = 0.55, main = "Module-trait correlation (FDR)")
dev.off()

## 6. Gene table and overlap with DEGs
read_de <- function(f) { d <- read.csv(f); rownames(d) <- d$gene_id; d }
is_de   <- function(d) !is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) >= 1
paired   <- read_de("results/DE_paired_tumor_vs_normal.csv")
unpaired <- read_de("results/DE_unpaired_tumor_vs_normal.csv")
de_both  <- intersect(paired$gene_id[is_de(paired)], unpaired$gene_id[is_de(unpaired)])

kme <- signedKME(dat, MEs)
own <- mapply(function(g, m) if (paste0("kME", m) %in% colnames(kme)) kme[g, paste0("kME", m)] else NA,
              colnames(dat), colors[colnames(dat)])
gt <- data.frame(gene_id = colnames(dat),
                 symbol  = gi$gene_name[match(colnames(dat), gi$gene_id)],
                 module  = colors[colnames(dat)],
                 kME     = as.numeric(own),
                 DE_both = colnames(dat) %in% de_both,
                 log2FC_paired = paired[colnames(dat), "log2FoldChange"],
                 padj_paired   = paired[colnames(dat), "padj"])
write.csv(gt[order(gt$module, -gt$kME), ], file.path(resdir, "wgcna_gene_modules.csv"), row.names = FALSE)

ov <- do.call(rbind, lapply(sort(unique(gt$module)), function(m) {
  a <- sum(gt$module == m & gt$DE_both);  b <- sum(gt$module == m & !gt$DE_both)
  c <- sum(gt$module != m & gt$DE_both);  d <- sum(gt$module != m & !gt$DE_both)
  data.frame(module = m, genes = a + b, DE_genes = a, DE_fraction = round(a / (a + b), 3),
             p = fisher.test(matrix(c(a, b, c, d), 2), alternative = "greater")$p.value)
}))
ov$FDR <- p.adjust(ov$p, "BH")
write.csv(ov[order(ov$FDR), ], file.path(resdir, "wgcna_module_deg_overlap.csv"), row.names = FALSE)

cat("\nDEGs (both analyses) among WGCNA genes:", sum(gt$DE_both), "of", nrow(gt), "\n")
cat("\nModule overlap with DEGs:\n"); print(ov[order(ov$FDR), ], row.names = FALSE)
cat("\nStrongest module-trait associations:\n")
print(head(long[order(long$FDR), ], 15), row.names = FALSE)
cat("\nSaved: results/wgcna_*.csv, figures/wgcna_*.png\n")
