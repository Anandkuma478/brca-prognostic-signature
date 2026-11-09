# 10_string_hubs.R
# Phase 2: STRING networks and hub genes for the DEG-enriched WGCNA modules.
# Hub rule (fixed before looking at any survival result): genes that are in the module and
# differentially expressed in both analyses; STRING v12 combined score >= 700;
# top 10 per module by Maximal Clique Centrality (MCC), ties broken by degree.
# Run from the project root with the brca_modules environment:
#   Rscript scripts/10_string_hubs.R 2>&1 | tee string_hubs.log
# Optional argument: output directory (default ".").

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."
resdir <- file.path(outdir, "results"); figdir <- file.path(outdir, "figures")
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
dir.create("data/string", showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(STRINGdb); library(igraph) })
set.seed(2031)
MODULES <- c("red", "green")
TOP_N   <- 10

genes <- read.csv("results/wgcna_gene_modules.csv")
sdb   <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 700,
                      network_type = "full", input_directory = "data/string")

log_mcc <- function(gr) {
  # log of MCC(v) = sum over maximal cliques containing v of (|C| - 1)!; logs avoid overflow
  cl   <- max_cliques(gr, min = 2)
  out  <- rep(-Inf, vcount(gr))
  size <- lengths(cl)
  for (v in seq_len(vcount(gr))) {
    has <- vapply(cl, function(x) v %in% as.integer(x), logical(1))
    if (any(has)) {
      lf <- lfactorial(size[has] - 1)
      out[v] <- max(lf) + log(sum(exp(lf - max(lf))))
    }
  }
  out / log(10)
}

all_nodes <- list(); all_edges <- list()
for (m in MODULES) {
  x <- genes[genes$module == m & genes$DE_both, ]
  mapped <- suppressWarnings(sdb$map(x, "symbol", removeUnmappedRows = TRUE))
  mapped <- mapped[!duplicated(mapped$STRING_id), ]
  e <- unique(sdb$get_interactions(mapped$STRING_id)[, c("from", "to", "combined_score")])
  sym <- setNames(mapped$symbol, mapped$STRING_id)
  e$from <- sym[e$from]; e$to <- sym[e$to]
  e <- e[!is.na(e$from) & !is.na(e$to) & e$from != e$to, ]

  gr <- simplify(graph_from_data_frame(e[, c("from", "to")], directed = FALSE,
                                       vertices = data.frame(name = mapped$symbol)))
  n <- data.frame(module = m, symbol = V(gr)$name, degree = degree(gr),
                  betweenness = round(betweenness(gr, normalized = TRUE), 4), log10_MCC = round(log_mcc(gr), 3))
  i <- match(n$symbol, x$symbol)
  n$kME <- round(x$kME[i], 3); n$log2FC_paired <- round(x$log2FC_paired[i], 3); n$padj_paired <- x$padj_paired[i]
  n <- n[order(-n$log10_MCC, -n$degree), ]
  n$hub <- seq_len(nrow(n)) <= TOP_N & n$degree > 0
  cat(sprintf("\n%s: %d DEGs, %d mapped to STRING, %d edges (score >= 700), %d genes with no edge\n",
              m, nrow(x), nrow(mapped), ecount(gr), sum(n$degree == 0)))
  print(head(n[, c("symbol", "degree", "log10_MCC", "kME", "log2FC_paired")], TOP_N), row.names = FALSE)
  all_nodes[[m]] <- n
  all_edges[[m]] <- cbind(module = m, e)

  hubs <- n$symbol[n$hub]
  sub  <- induced_subgraph(gr, hubs)
  png(file.path(figdir, sprintf("string_hubs_%s.png", m)), width = 1500, height = 1500, res = 220)
  par(mar = c(1, 1, 3, 1))
  plot(sub, layout = layout_in_circle(sub), vertex.size = 30, vertex.color = m, vertex.frame.color = "grey30",
       vertex.label.color = if (m == "red") "white" else "black", vertex.label.cex = 0.8,
       edge.color = "grey60", main = sprintf("Top %d hub genes, %s module", length(hubs), m))
  dev.off()
}
nodes <- do.call(rbind, all_nodes); edges <- do.call(rbind, all_edges)
rownames(nodes) <- NULL
write.csv(nodes, file.path(resdir, "string_hub_genes.csv"), row.names = FALSE)
write.csv(edges, file.path(resdir, "string_edges.csv"), row.names = FALSE)
cat("\nSaved: results/string_hub_genes.csv, results/string_edges.csv, figures/string_hubs_*.png\n")
