# 11_enrichment.R
# Phase 2: GO and KEGG over-representation for each WGCNA module, and Hallmark GSEA on the paired DE results.
# Background for over-representation: the genes entered into WGCNA.
# Run from the project root with the brca_modules environment:
#   Rscript scripts/11_enrichment.R 2>&1 | tee enrichment.log
# Optional argument: output directory (default ".").

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."
resdir <- file.path(outdir, "results"); figdir <- file.path(outdir, "figures")
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr); library(ggplot2)
})
set.seed(2032)

genes <- read.csv("results/wgcna_gene_modules.csv")
map   <- suppressWarnings(bitr(genes$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
map   <- map[!duplicated(map$SYMBOL), ]
genes$entrez <- map$ENTREZID[match(genes$symbol, map$SYMBOL)]
universe <- unique(na.omit(genes$entrez))
cat("WGCNA genes:", nrow(genes), " Mapped to Entrez (background):", length(universe), "\n")

modules <- setdiff(sort(unique(genes$module)), "grey")
go_all <- list(); kegg_all <- list(); kegg_failed <- FALSE
for (m in modules) {
  ids <- unique(na.omit(genes$entrez[genes$module == m]))
  ego <- enrichGO(gene = ids, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
    d <- as.data.frame(ego)
    if (nrow(d)) go_all[[m]] <- cbind(module = m, d)
  }
  ek <- tryCatch(enrichKEGG(gene = ids, universe = universe, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05),
                 error = function(e) { kegg_failed <<- TRUE; message("KEGG failed for ", m, ": ", conditionMessage(e)); NULL })
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) kegg_all[[m]] <- cbind(module = m, as.data.frame(setReadable(ek, org.Hs.eg.db, "ENTREZID")))
  cat(sprintf("%s: %d genes | GO terms %d | KEGG pathways %d\n", m, length(ids),
              if (is.null(go_all[[m]])) 0L else nrow(go_all[[m]]), if (is.null(kegg_all[[m]])) 0L else nrow(kegg_all[[m]])))
}
go <- do.call(rbind, go_all); rownames(go) <- NULL
write.csv(go, file.path(resdir, "enrichment_go_bp.csv"), row.names = FALSE)
if (length(kegg_all)) {
  kegg <- do.call(rbind, kegg_all); rownames(kegg) <- NULL
  write.csv(kegg, file.path(resdir, "enrichment_kegg.csv"), row.names = FALSE)
}
if (kegg_failed) cat("WARNING: KEGG could not be reached for at least one module.\n")

top_terms <- function(d, n) do.call(rbind, lapply(split(d, d$module), function(x) head(x[order(x$p.adjust), ], n)))
cat("\nTop GO BP terms per module:\n")
print(top_terms(go, 5)[, c("module", "Description", "Count", "p.adjust")], row.names = FALSE)
if (length(kegg_all)) { cat("\nTop KEGG pathways per module:\n"); print(top_terms(kegg, 5)[, c("module", "Description", "Count", "p.adjust")], row.names = FALSE) }

p <- top_terms(go, 5)
p$Description <- substr(p$Description, 1, 55)
p$label <- factor(paste0(p$Description, " [", p$module, "]"), levels = rev(unique(paste0(p$Description, " [", p$module, "]"))))
g <- ggplot(p, aes(-log10(p.adjust), label, size = Count, colour = module)) + geom_point() +
  scale_colour_manual(values = setNames(modules, modules), guide = "none") +
  labs(x = "-log10 adjusted p", y = NULL, size = "Genes", title = "GO biological process, by WGCNA module") +
  theme_bw(base_size = 10)
ggsave(file.path(figdir, "enrichment_go_modules.png"), g, width = 9, height = 7, dpi = 600)

## Hallmark GSEA, ranked by the paired Wald statistic
de <- read.csv("results/DE_paired_tumor_vs_normal.csv")
de <- de[!is.na(de$stat) & !is.na(de$gene_name) & de$gene_name != "", ]
de <- de[order(-abs(de$stat)), ]; de <- de[!duplicated(de$gene_name), ]
ranks <- sort(setNames(de$stat, de$gene_name), decreasing = TRUE)
h <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "H"))[, c("gs_name", "gene_symbol")]
gs <- GSEA(ranks, TERM2GENE = h, minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)
gd <- as.data.frame(gs)[, c("ID", "setSize", "NES", "pvalue", "p.adjust")]
gd <- gd[order(gd$p.adjust), ]; rownames(gd) <- NULL
write.csv(gd, file.path(resdir, "enrichment_gsea_hallmark.csv"), row.names = FALSE)
cat("\nHallmark GSEA (paired tumour vs normal): sets tested", nrow(gd), " FDR<0.05:", sum(gd$p.adjust < 0.05), "\n")
print(head(gd, 15), row.names = FALSE)

s <- head(gd[gd$p.adjust < 0.05, ], 20)
s$ID <- factor(sub("HALLMARK_", "", s$ID), levels = sub("HALLMARK_", "", s$ID[order(s$NES)]))
g2 <- ggplot(s, aes(NES, ID, fill = NES > 0)) + geom_col() +
  scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
  labs(x = "Normalised enrichment score (tumour vs normal)", y = NULL, title = "Hallmark GSEA") + theme_bw(base_size = 10)
ggsave(file.path(figdir, "enrichment_gsea_hallmark.png"), g2, width = 8, height = 6, dpi = 600)
cat("\nSaved: results/enrichment_*.csv, figures/enrichment_*.png\n")
