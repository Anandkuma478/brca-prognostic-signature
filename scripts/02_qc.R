library(SummarizedExperiment)
library(edgeR)
library(DESeq2)
library(ggplot2)

setwd("/home/administrator/Desktop/anand 2/brca genomics")

se   <- readRDS("data/tcga_brca_se.rds")
meta <- as.data.frame(colData(se))
cat("Sample types:\n"); print(table(meta$sample_type))

keep_s <- meta$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
se   <- se[, keep_s]
meta <- as.data.frame(colData(se))

meta$patient <- substr(meta$barcode, 1, 12)
# Count repeats within sample type only: a patient's normal must not mark their tumor as a duplicate
is_dup <- meta$sample_type == "Primary Tumor" & duplicated(paste(meta$patient, meta$sample_type))
se   <- se[, !is_dup]
meta <- as.data.frame(colData(se))

meta$group <- factor(ifelse(meta$sample_type == "Primary Tumor", "Tumor", "Normal"),
                     levels = c("Normal", "Tumor"))
counts <- assay(se, "unstranded")

keep <- filterByExpr(counts, group = meta$group)
counts <- counts[keep, ]
cat("After filtering: genes =", nrow(counts), " samples =", ncol(counts), "\n")
print(table(meta$group))

gene_info <- as.data.frame(rowData(se))[rownames(counts), c("gene_id", "gene_name", "gene_type")]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ group)
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))

p <- ggplot(pca, aes(PC1, PC2, color = group)) +
  geom_point(alpha = 0.7) +
  xlab(paste0("PC1: ", pv[1], "%")) + ylab(paste0("PC2: ", pv[2], "%")) +
  theme_bw()
ggsave("figures/pca_tumor_vs_normal.png", p, width = 6, height = 5, dpi = 300)

saveRDS(counts,    "data/counts_filtered.rds")
saveRDS(meta,      "data/meta_clean.rds")
saveRDS(gene_info, "data/gene_info.rds")
cat("Done. PCA plot: figures/pca_tumor_vs_normal.png\n")
