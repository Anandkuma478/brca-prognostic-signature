# 03b_volcano_heatmap.R
# Volcano plots (paired + unpaired DE) and a tumor-vs-normal heatmap of the top DEGs.
# Thresholds from the plan: |log2FC| >= 1 and adjusted p < 0.05.
# Run from the project root:
#   Rscript scripts/03b_volcano_heatmap.R 2>&1 | tee volcano_heatmap.log
# Outputs: figures/volcano_paired.(pdf|png), figures/volcano_unpaired.(pdf|png),
#          figures/heatmap_top50_paired.pdf, results/DEG_counts.csv

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
})
dir.create("figures", showWarnings = FALSE)

LFC  <- 1
PADJ <- 0.05
TOPN <- 25                       # per direction in the heatmap (25 up + 25 down)
HAS_REPEL <- requireNamespace("ggrepel", quietly = TRUE)

pick <- function(nms, cands) { h <- intersect(cands, nms); if (length(h)) h[1] else NA_character_ }

# ---------------------------------------------------------------- load DE
load_de <- function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  cat("\n", f, "\n  columns:", paste(names(d), collapse = ", "), "\n")
  lfc <- pick(names(d), c("log2FoldChange", "logFC", "log2FC"))
  p   <- pick(names(d), c("padj", "adj.P.Val", "FDR", "adj_p"))
  gid <- pick(names(d), c("gene_id", "ensembl_id", "Geneid", "X", names(d)[1]))
  sym <- pick(names(d), c("gene_name", "symbol", "SYMBOL", "gene_symbol"))
  if (is.na(lfc) || is.na(p))
    stop("Could not find log2FC / adjusted-p columns in ", f, ". Columns are listed above; send them to me.")
  out <- data.frame(gene_id = as.character(d[[gid]]),
                    symbol  = if (!is.na(sym)) as.character(d[[sym]]) else as.character(d[[gid]]),
                    lfc = d[[lfc]], padj = d[[p]], stringsAsFactors = FALSE)
  # fall back to data/gene_info.rds for symbols if the DE file has none
  if (is.na(sym)) {
    gi <- tryCatch(as.data.frame(readRDS("data/gene_info.rds")), error = function(e) NULL)
    if (!is.null(gi)) {
      idc <- pick(names(gi), c("gene_id", "ensembl_gene_id"))
      nmc <- pick(names(gi), c("gene_name", "symbol", "external_gene_name"))
      if (!is.na(idc) && !is.na(nmc)) {
        m <- match(out$gene_id, gi[[idc]])
        out$symbol <- ifelse(is.na(m) | gi[[nmc]][m] == "", out$gene_id, gi[[nmc]][m])
      }
    }
  }
  out[!is.na(out$lfc) & !is.na(out$padj), ]
}

classify <- function(d) {
  ifelse(d$padj < PADJ & d$lfc >= LFC, "Up",
         ifelse(d$padj < PADJ & d$lfc <= -LFC, "Down", "NS"))
}

# ---------------------------------------------------------------- volcano
volcano <- function(d, title, stem) {
  d$padj[d$padj < 1e-300] <- 1e-300
  d$class <- classify(d)
  n <- table(factor(d$class, levels = c("Up", "Down", "NS")))
  cat(sprintf("%s: Up = %d, Down = %d, NS = %d (of %d genes)\n", title, n["Up"], n["Down"], n["NS"], nrow(d)))
  sig <- d[d$class != "NS", ]
  top <- head(sig[order(sig$padj), ], 15)
  g <- ggplot(d, aes(lfc, -log10(padj), colour = class)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_colour_manual(values = c(Down = "#2166ac", NS = "grey75", Up = "#b2182b"),
                        breaks = c("Up", "Down", "NS")) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(PADJ), linetype = "dashed", colour = "grey40") +
    labs(title = title,
         subtitle = sprintf("Up %d | Down %d  (|log2FC| >= %g, adj. p < %g)", n["Up"], n["Down"], LFC, PADJ),
         x = "log2 fold change (tumor vs normal)", y = "-log10 adjusted p", colour = NULL) +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
  if (nrow(top)) {
    g <- g + if (HAS_REPEL)
      ggrepel::geom_text_repel(data = top, aes(label = symbol), size = 3, colour = "black",
                               max.overlaps = Inf, show.legend = FALSE)
    else geom_text(data = top, aes(label = symbol), size = 2.5, colour = "black",
                   vjust = -0.5, show.legend = FALSE)
  }
  ggsave(paste0("figures/", stem, ".pdf"), g, width = 7, height = 6)
  ggsave(paste0("figures/", stem, ".png"), g, width = 7, height = 6, dpi = 300)
  data.frame(analysis = title, up = n[["Up"]], down = n[["Down"]], ns = n[["NS"]], total = nrow(d))
}

paired   <- load_de("results/DE_paired_tumor_vs_normal.csv")
unpaired <- load_de("results/DE_unpaired_tumor_vs_normal.csv")
m0 <- as.data.frame(readRDS("data/meta_clean.rds"))
n_pairs <- length(intersect(m0$patient[m0$shortLetterCode == "TP"], m0$patient[m0$shortLetterCode == "NT"]))
counts_tbl <- rbind(volcano(paired,   sprintf("Paired tumor vs normal (%d pairs)", n_pairs), "volcano_paired"),
                    volcano(unpaired, "Unpaired tumor vs normal",          "volcano_unpaired"))
write.csv(counts_tbl, "results/DEG_counts.csv", row.names = FALSE)

# ---------------------------------------------------------------- heatmap (paired samples)
meta <- as.data.frame(readRDS("data/meta_clean.rds"))
rownames(meta) <- meta$barcode
cnt <- readRDS("data/counts_filtered.rds")
if (inherits(cnt, "SummarizedExperiment")) cnt <- SummarizedExperiment::assay(cnt)
cnt <- as.matrix(cnt)

meta <- meta[meta$barcode %in% colnames(cnt), ]
tp <- meta[meta$shortLetterCode == "TP", ]
nt <- meta[meta$shortLetterCode == "NT", ]
pp <- intersect(tp$patient, nt$patient)
tp <- tp[!duplicated(tp$patient) & tp$patient %in% pp, ]
nt <- nt[!duplicated(nt$patient) & nt$patient %in% pp, ]
cat(sprintf("\nHeatmap samples: %d tumor + %d normal (paired patients: %d)\n", nrow(tp), nrow(nt), length(pp)))

samp <- c(tp$barcode, nt$barcode)
grp  <- factor(c(rep("Tumor", nrow(tp)), rep("Normal", nrow(nt))), levels = c("Normal", "Tumor"))
x <- cnt[, samp, drop = FALSE]
if (max(x, na.rm = TRUE) > 1000) {               # raw counts -> log2 CPM
  x <- log2(t(t(x) / colSums(x) * 1e6) + 1)
}

sig <- paired[paired$padj < PADJ & abs(paired$lfc) >= LFC & paired$gene_id %in% rownames(x), ]
if (!nrow(sig)) stop("No paired DEGs match the rownames of counts_filtered.rds; gene ID format differs. Send me rownames(counts)[1:3] and head(DE).")
up   <- head(sig[sig$lfc > 0, ][order(sig$padj[sig$lfc > 0]), ], TOPN)
down <- head(sig[sig$lfc < 0, ][order(sig$padj[sig$lfc < 0]), ], TOPN)
hg <- rbind(up, down)

m <- x[hg$gene_id, , drop = FALSE]
m <- t(scale(t(m)))                               # z-score per gene
m[m >  3] <-  3; m[m < -3] <- -3
rownames(m) <- make.unique(hg$symbol)
m <- m[, order(grp), drop = FALSE]
ann <- data.frame(Group = grp[order(grp)], row.names = colnames(m))

pheatmap(m, annotation_col = ann, cluster_cols = FALSE, cluster_rows = TRUE,
         show_colnames = FALSE, fontsize_row = 7,
         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         annotation_colors = list(Group = c(Normal = "#2166ac", Tumor = "#b2182b")),
         main = sprintf("Top %d up + %d down DEGs (paired, %d pairs)", nrow(up), nrow(down), length(pp)),
         filename = "figures/heatmap_top50_paired.pdf", width = 8, height = 9)

cat("\nSaved: figures/volcano_paired.*, figures/volcano_unpaired.*, figures/heatmap_top50_paired.pdf, results/DEG_counts.csv\n")
print(counts_tbl)
