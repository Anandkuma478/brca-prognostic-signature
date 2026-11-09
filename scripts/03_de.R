suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

counts <- readRDS("data/counts_filtered.rds")
if (inherits(counts, "SummarizedExperiment")) counts <- assay(counts, "unstranded")
counts <- as.matrix(counts)
meta   <- as.data.frame(readRDS("data/meta_clean.rds"))
ginfo  <- tryCatch(as.data.frame(readRDS("data/gene_info.rds")), error = function(e) NULL)

if (!all(colnames(counts) %in% rownames(meta))) stop("colnames(counts) != rownames(meta)")
meta <- meta[colnames(counts), , drop = FALSE]

if (any(counts != round(counts))) message("Non-integer values: rounding. Check these are raw counts.")
counts <- round(counts); mode(counts) <- "integer"

bc <- colnames(counts)
meta$patient   <- substr(bc, 1, 12)
meta$condition <- factor(ifelse(substr(bc, 14, 15) == "11", "Normal", "Tumor"),
                         levels = c("Normal", "Tumor"))
cat("Samples:", table(meta$condition), "\n")

annotate <- function(df) {
  df$gene_id <- rownames(df)
  if (!is.null(ginfo)) {
    id_col <- intersect(c("gene_id", "ensembl_gene_id", "id"), colnames(ginfo))[1]
    nm_col <- intersect(c("gene_name", "external_gene_name", "symbol"), colnames(ginfo))[1]
    if (!is.na(id_col) && !is.na(nm_col)) {
      df$gene_name <- ginfo[[nm_col]][match(df$gene_id, ginfo[[id_col]])]
    } else if (!is.na(nm_col)) {
      df$gene_name <- ginfo[[nm_col]][match(df$gene_id, rownames(ginfo))]
    }
  }
  df[order(df$padj), ]
}

volcano <- function(df, file, title) {
  df <- df[!is.na(df$padj), ]
  df$sig <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) >= 1,
                   ifelse(df$log2FoldChange > 0, "Up in Tumor", "Down in Tumor"), "NS")
  cols <- c("Up in Tumor" = "#d62728", "Down in Tumor" = "#1f77b4", "NS" = "grey70")
  png(file, width = 1800, height = 1500, res = 250)
  plot(df$log2FoldChange, -log10(df$padj), pch = 16, cex = 0.4, col = cols[df$sig],
       xlab = "log2 fold change (Tumor vs Normal)", ylab = "-log10 adjusted p", main = title)
  abline(v = c(-1, 1), h = -log10(0.05), lty = 2, col = "grey40")
  legend("topright", legend = names(cols), col = cols, pch = 16, bty = "n")
  dev.off()
}

## Paired (patients with both tumor and normal)
paired_pts <- intersect(meta$patient[meta$condition == "Tumor"],
                        meta$patient[meta$condition == "Normal"])
cat("Paired patients:", length(paired_pts), "\n")
keep <- meta$patient %in% paired_pts
mp <- meta[keep, ]; mp$patient <- factor(mp$patient)
dds_p <- DESeqDataSetFromMatrix(counts[, keep], colData = mp, design = ~ patient + condition)
dds_p <- DESeq(dds_p)
res_p <- tryCatch(lfcShrink(dds_p, coef = "condition_Tumor_vs_Normal", type = "apeglm"),
                  error = function(e) { message("apeglm not available, using unshrunk results"); results(dds_p, name = "condition_Tumor_vs_Normal") })
res_p <- annotate(as.data.frame(res_p))
write.csv(res_p, "results/DE_paired_tumor_vs_normal.csv", row.names = FALSE)
volcano(res_p, "figures/volcano_paired.png", sprintf("Paired: Tumor vs Normal (%d pairs)", length(paired_pts)))

## Unpaired (all samples)
dds_u <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ condition)
dds_u <- DESeq(dds_u)
res_u <- annotate(as.data.frame(results(dds_u, name = "condition_Tumor_vs_Normal")))
write.csv(res_u, "results/DE_unpaired_tumor_vs_normal.csv", row.names = FALSE)
volcano(res_u, "figures/volcano_unpaired.png", "All samples: Tumor vs Normal")

## Summary + BRCA1/2
summ <- function(df, label) {
  s <- df[!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) >= 1, ]
  cat(label, ": ", nrow(s), " DE genes; up=", sum(s$log2FoldChange > 0),
      " down=", sum(s$log2FoldChange < 0), "\n", sep = "")
}
summ(res_p, "Paired"); summ(res_u, "Unpaired")
if ("gene_name" %in% colnames(res_p)) {
  cols <- c("gene_id", "gene_name", "baseMean", "log2FoldChange", "padj")
  cat("\nBRCA1/BRCA2 paired:\n");   print(res_p[res_p$gene_name %in% c("BRCA1", "BRCA2"), cols])
  cat("\nBRCA1/BRCA2 unpaired:\n"); print(res_u[res_u$gene_name %in% c("BRCA1", "BRCA2"), cols])
}
cat("\nDone. Results in results/, plots in figures/\n")
