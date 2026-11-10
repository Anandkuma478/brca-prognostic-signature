# 21_final_figures.R
# Draws manuscript figures at 600 dpi from saved result tables. No analysis is repeated, so every
# number shown is exactly the one in the result table.
#   Figure 5:  module-trait correlation heatmap (results/wgcna_module_trait.csv, from 08_wgcna.R)
#   Figure 6:  STRING networks of the hub genes  (results/string_hub_genes.csv, string_edges.csv, from 10_string_hubs.R)
#   Figure 7:  GO dot plot                       (results/enrichment_go_bp.csv, from 11_enrichment.R)
#   Figure 8:  Hallmark GSEA bar plot             (results/enrichment_gsea_hallmark.csv, from 11_enrichment.R)
#   Figure 12: immune heatmap                     (results/immune_association.csv, from 16_immune_association.R)
#   Figure 13: mutation dot plot                  (results/mutation_association.csv, from 18_mutation_association.R)
#   Figure 14: published-signature forest plot    (results/published_signature_comparison.csv, from 20_published_and_random_benchmark.R)
# Figures 7, 8, 12, 13 and 14 use the same plotting code as the scripts named above, which run online
# queries (KEGG) or long resampling, so they are redrawn here from the saved tables rather than by rerunning.
# Run from the project root with the brca_modules environment:
#   Rscript scripts/21_final_figures.R

suppressPackageStartupMessages({ library(ggplot2); library(igraph) })
dir.create("figures", showWarnings = FALSE)
DPI <- 600

## Figure 5: module-trait correlations
mt <- read.csv("results/wgcna_module_trait.csv")
mod_order <- c("green", "red", "brown", "turquoise", "blue", "yellow", "grey")
mod_lab <- setNames(c("Green", "Red", "Brown", "Turquoise", "Blue", "Yellow", "Grey (unassigned)"), mod_order)
trait_lab <- c(Age = "Age", Stage = "Stage", LumA = "Luminal A", LumB = "Luminal B", Her2 = "HER2-enriched",
               Basal = "Basal-like", Normal_like = "Normal-like")
stopifnot(setequal(mt$module, mod_order), setequal(mt$trait, names(trait_lab)))
fmt_fdr <- function(q) ifelse(q < 0.001, sub("e-0?", "e-", sprintf("%.0e", q)), sprintf("%.3f", q))
mt$Module <- factor(mod_lab[mt$module], levels = rev(mod_lab))
mt$Trait  <- factor(trait_lab[mt$trait], levels = trait_lab)
mt$label  <- sprintf("%.2f\n(%s)", mt$cor, fmt_fdr(mt$FDR))
g <- ggplot(mt, aes(Trait, Module, fill = cor)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label, colour = abs(cor) > 0.55), size = 3.1, lineheight = 0.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#222222")) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, limits = c(-1, 1),
                       name = "Pearson r") +
  labs(x = NULL, y = "Module eigengene") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1, colour = "#222222"),
        axis.text.y = element_text(colour = "#222222"), legend.title = element_text(size = 10))
ggsave("figures/fig_module_trait.png", g, width = 7.2, height = 5.2, dpi = DPI, bg = "white")

## Figure 6: STRING interactions among the hub genes
nodes <- read.csv("results/string_hub_genes.csv")
edges <- read.csv("results/string_edges.csv")
fill <- c(red = "#f3b3b3", green = "#b5dfb5")
png("figures/fig_string_hubs.png", width = 7.2, height = 3.7, units = "in", res = DPI)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.2, 0.5))
for (m in c("red", "green")) {
  hubs <- nodes$symbol[nodes$module == m & nodes$hub]
  e <- edges[edges$module == m & edges$from %in% hubs & edges$to %in% hubs, c("from", "to")]
  gr <- graph_from_data_frame(e, directed = FALSE, vertices = data.frame(name = hubs))
  plot(gr, layout = layout_in_circle(gr), vertex.size = 42, vertex.color = fill[[m]], vertex.frame.color = "grey35",
       vertex.label.color = "black", vertex.label.cex = 0.6, vertex.label.family = "sans", vertex.label.font = 2,
       edge.color = "grey55", edge.width = 1.1, margin = 0.08)
  title(sprintf("%s module (%d interactions)", tools::toTitleCase(m), nrow(e)), cex.main = 1.05, font.main = 1)
}
invisible(dev.off())

## Figure 7: GO biological process terms per module (plotting code as in 11_enrichment.R)
genes <- read.csv("results/wgcna_gene_modules.csv")
modules <- setdiff(sort(unique(genes$module)), "grey")
go <- read.csv("results/enrichment_go_bp.csv")
top_terms <- function(d, n) do.call(rbind, lapply(split(d, d$module), function(x) head(x[order(x$p.adjust), ], n)))
p <- top_terms(go, 5)
p$Description <- substr(p$Description, 1, 55)
p$label <- factor(paste0(p$Description, " [", p$module, "]"), levels = rev(unique(paste0(p$Description, " [", p$module, "]"))))
g <- ggplot(p, aes(-log10(p.adjust), label, size = Count, colour = module)) + geom_point() +
  scale_colour_manual(values = setNames(modules, modules), guide = "none") +
  labs(x = "-log10 adjusted p", y = NULL, size = "Genes", title = "GO biological process, by WGCNA module") +
  theme_bw(base_size = 10)
ggsave("figures/enrichment_go_modules.png", g, width = 9, height = 7, dpi = DPI)

## Figure 8: Hallmark GSEA (plotting code as in 11_enrichment.R)
gd <- read.csv("results/enrichment_gsea_hallmark.csv")
gd <- gd[order(gd$p.adjust), ]
s <- head(gd[gd$p.adjust < 0.05, ], 20)
s$ID <- factor(sub("HALLMARK_", "", s$ID), levels = sub("HALLMARK_", "", s$ID[order(s$NES)]))
g2 <- ggplot(s, aes(NES, ID, fill = NES > 0)) + geom_col() +
  scale_fill_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
  labs(x = "Normalised enrichment score (tumour vs normal)", y = NULL, title = "Hallmark GSEA") + theme_bw(base_size = 10)
ggsave("figures/enrichment_gsea_hallmark.png", g2, width = 8, height = 6, dpi = DPI)

## Figure 12: immune features (plotting code as in 16_immune_association.R)
SETS <- paste0("HALLMARK_", c("INTERFERON_GAMMA_RESPONSE", "INTERFERON_ALPHA_RESPONSE", "INFLAMMATORY_RESPONSE", "ALLOGRAFT_REJECTION",
                              "IL2_STAT5_SIGNALING", "IL6_JAK_STAT3_SIGNALING", "TNFA_SIGNALING_VIA_NFKB", "COMPLEMENT"))
CHECK <- c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1")
res <- read.csv("results/immune_association.csv")
pretty <- function(x) { x <- gsub("_", " ", x); x <- sub("Cytolytic score", "Cytolytic score (GZMA, PRF1)", x); x }
order_f <- c("Cytolytic_score", "CD8A", CHECK, sub("HALLMARK_", "", SETS))
stopifnot(all(res$feature %in% order_f))
res$feature_lab <- factor(pretty(res$feature), levels = rev(pretty(order_f)))
res$cohort <- factor(res$cohort, levels = c("SCAN-B", "TCGA-BRCA", "METABRIC"))
res$panel <- paste("Model", res$model)
g <- ggplot(res, aes(cohort, feature_lab, fill = rho)) + geom_tile(colour = "#fcfcfb", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.4, colour = "#0b0b0b") + facet_wrap(~ panel) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-0.8, 0.8), oob = scales::squish, name = "Spearman rho") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"), axis.text.y = element_text(colour = "#0b0b0b"))
ggsave("figures/immune_association_heatmap.png", g, width = 8, height = 5.6, dpi = DPI, bg = "white")

## Figure 13: mutation frequency by risk group (plotting code as in 18_mutation_association.R)
res <- read.csv("results/mutation_association.csv")
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
ggsave("figures/mutation_by_risk.png", g, width = 7.5, height = 5 + 0.08 * length(common), dpi = DPI, bg = "white")

## Figure 14: C-index of the models and published signatures (plotting code as in 20_published_and_random_benchmark.R)
res <- read.csv("results/published_signature_comparison.csv")
f <- res[!is.na(res$c_index), ]
f$cohort <- factor(f$cohort, levels = c("SCAN-B (out-of-fold)", "TCGA-BRCA", "METABRIC"))
f$signature <- factor(f$signature, levels = rev(c("Model A", "Model B", "Oncotype DX (21-gene)", "MammaPrint (70-gene)", "PAM50 ROR-S", "Genomic Grade Index")))
f$Group <- ifelse(grepl("^Model", f$signature), "This study", "Published (research-based)")
g <- ggplot(f, aes(c_index, signature, colour = Group)) + geom_vline(xintercept = 0.5, linetype = 2, colour = "#c3c2b7") +
  geom_segment(aes(x = c_low, xend = c_high, yend = signature), linewidth = 0.7) + geom_point(size = 2.2) +
  geom_text(aes(x = 0.81, label = sprintf("%.3f", c_index)), colour = "#52514e", size = 3.3, hjust = 1, show.legend = FALSE) +
  facet_wrap(~ cohort) + scale_colour_manual(values = c(`This study` = "#2a78d6", `Published (research-based)` = "#eb6834"), name = NULL) +
  scale_x_continuous(limits = c(0.45, 0.82), breaks = c(0.5, 0.6, 0.7, 0.8)) + labs(x = "C-index for overall survival (95% bootstrap CI)", y = NULL) + theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = "#ecebe7", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "top", panel.border = element_rect(colour = "#c3c2b7"))
ggsave("figures/published_signature_forest.png", g, width = 8, height = 3.3, dpi = DPI, bg = "white")

cat("Saved at", DPI, "dpi: figures/fig_module_trait.png, fig_string_hubs.png, enrichment_go_modules.png, enrichment_gsea_hallmark.png,",
    "immune_association_heatmap.png, mutation_by_risk.png, published_signature_forest.png\n")
