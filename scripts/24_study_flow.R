# 24_study_flow.R
# Study flow diagram (REMARK item 12): patients and genes at each stage of the analysis.
# Counts are those reported by 01-06c and 09 (see the corresponding logs).
# Run from the project root with the brca_modules environment:
#   Rscript scripts/24_study_flow.R

suppressPackageStartupMessages(library(ggplot2))

box <- function(id, x, y, w, h, label, fill = "#f3f2ee") data.frame(id, x, y, w, h, label, fill)
b <- rbind(
  box("t1", 2.55, 11.0, 4.9, 1.25, "TCGA-BRCA (Genomic Data Commons)\n1231 samples: 1111 primary tumour,\n113 solid tissue normal, 7 metastatic"),
  box("t2", 2.55, 9.15, 4.9, 1.25, "Quality control\nExcluded: 7 metastatic, 16 repeat tumour samples\n1095 primary tumours, 113 normal samples"),
  box("t3", 2.55, 7.3, 4.9, 1.25, "Tumour vs normal differential expression\nPaired (113 pairs) and unpaired (1208 samples)\n5068 genes significant in both analyses"),
  box("t4", 2.55, 5.45, 4.9, 1.25, "Candidate filters (no survival data used)\nProtein-coding, present and expressed in SCAN-B\n1688 candidate genes"),
  box("s1", 7.75, 7.3, 4.3, 1.25, "SCAN-B (GEO GSE96058)\n136 technical replicates removed\n3273 patients, 336 deaths"),
  box("m",  5.15, 3.6, 9.9, 1.25, "Model development in SCAN-B: LASSO Cox regression, nested cross-validation\nModel A (85 genes) and Model B (29 genes)\nCoefficients and risk-group cut-off frozen", "#e3edf8"),
  box("v1", 2.55, 1.0, 4.9, 1.45, "External validation: TCGA-BRCA\n1071 patients, 149 deaths\n24 excluded (missing or non-positive survival)\nSame tumours as used for candidate selection", "#e6f2ea"),
  box("v2", 7.75, 1.0, 4.9, 1.45, "External validation: METABRIC (cBioPortal)\n1979 patients, 1143 deaths\n1 excluded (missing or non-positive survival)\nNot used before validation", "#e6f2ea"))
top <- function(i) c(b$x[b$id == i], b$y[b$id == i] + b$h[b$id == i] / 2)
bot <- function(i) c(b$x[b$id == i], b$y[b$id == i] - b$h[b$id == i] / 2)
arrow_df <- function(from, to) data.frame(x = from[1], y = from[2], xend = to[1], yend = to[2])
mtop <- b$y[b$id == "m"] + b$h[b$id == "m"] / 2; mbot <- b$y[b$id == "m"] - b$h[b$id == "m"] / 2
a <- rbind(arrow_df(bot("t1"), top("t2")), arrow_df(bot("t2"), top("t3")), arrow_df(bot("t3"), top("t4")),
           arrow_df(bot("t4"), c(2.55, mtop)), arrow_df(bot("s1"), c(7.75, mtop)),
           arrow_df(c(2.55, mbot), top("v1")), arrow_df(c(7.75, mbot), top("v2")))

g <- ggplot() +
  geom_segment(data = a, aes(x, y, xend = xend, yend = yend), linewidth = 0.45, colour = "#444444",
               arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
  geom_rect(data = b, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = fill),
            colour = "#555555", linewidth = 0.4) +
  geom_text(data = b, aes(x, y, label = label), size = 2.85, lineheight = 0.95, family = "Arial") +
  scale_fill_identity() + coord_cartesian(xlim = c(-0.05, 10.25), ylim = c(0.1, 11.7), expand = FALSE) +
  theme_void()
ggsave("figures/fig_study_flow.png", g, width = 6.85, height = 6.6, dpi = 650, bg = "white")
cat("Saved: figures/fig_study_flow.png\n")
