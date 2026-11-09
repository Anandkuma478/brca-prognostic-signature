# 19_protein_atlas.R
# Phase 4: descriptive protein-level support from Human Protein Atlas immunohistochemistry.
# Genes: the ten signature genes selected in at least 50% of bootstrap resamples, and the 20 hub genes.
# Input files in data/hpa/: cancer_data.tsv and normal_ihc_data.tsv (www.proteinatlas.org, CC BY-SA 4.0).
# Run from the project root:  Rscript scripts/19_protein_atlas.R 2>&1 | tee protein.log

st <- read.csv("results/stability_selection.csv"); stable <- st$symbol[st$freq >= 0.5]
hubs <- read.csv("results/string_hub_genes.csv"); hubs <- hubs[hubs$hub, ]
cf <- read.csv("results/signature_scanb_sparse30.csv")
de <- read.csv("results/DE_paired_tumor_vs_normal.csv")
genes <- data.frame(symbol = c(stable, hubs$symbol),
                    group = c(rep("Stable signature gene", length(stable)), paste("Hub gene,", hubs$module, "module")))

ca <- read.delim("data/hpa/cancer_data.tsv", quote = "", check.names = FALSE)
ca <- ca[ca$Cancer == "breast cancer", ]
no <- read.delim("data/hpa/normal_ihc_data.tsv", quote = "", check.names = FALSE)
no <- no[no$Tissue == "Breast" & no$`Cell type` == "glandular cells", ]

out <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
  g <- genes$symbol[i]; a <- ca[ca$`Gene name` == g, ]; b <- no[no$`Gene name` == g, ]
  n <- if (nrow(a)) sum(a[1, c("High", "Medium", "Low", "Not detected")]) else 0L
  d <- de[which(de$gene_name == g)[1], ]
  data.frame(gene = g, group = genes$group[i],
             model_B_coefficient = if (g %in% cf$symbol) cf$coef[cf$symbol == g] else NA,
             paired_log2FC = d$log2FoldChange,
             tumours_stained = n,
             high = if (n) a$High[1] else NA, medium = if (n) a$Medium[1] else NA, low = if (n) a$Low[1] else NA,
             not_detected = if (n) a$`Not detected`[1] else NA,
             detected_pct = if (n) round(100 * (a$High[1] + a$Medium[1] + a$Low[1]) / n) else NA,
             normal_glandular_level = if (nrow(b)) b$Level[1] else NA,
             normal_reliability = if (nrow(b)) b$Reliability[1] else NA)
}))
write.csv(out, "results/protein_atlas_support.csv", row.names = FALSE)
print(out, row.names = FALSE)
cat("\nGenes with tumour staining data:", sum(out$tumours_stained > 0), "of", nrow(out), "\n")
cat("Saved: results/protein_atlas_support.csv\n")
