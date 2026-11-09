# 17_tcga_maf_download.R
# Download open-access masked somatic mutation calls for TCGA-BRCA from the Genomic Data Commons.
# Run from the project root with the brca_genomics environment:
#   Rscript scripts/17_tcga_maf_download.R 2>&1 | tee tcga_maf.log

suppressPackageStartupMessages(library(TCGAbiolinks))
dir.create("data/tcga_maf", showWarnings = FALSE, recursive = TRUE)
q <- GDCquery(project = "TCGA-BRCA", data.category = "Simple Nucleotide Variation", access = "open",
              data.type = "Masked Somatic Mutation", workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(q, directory = "data/tcga_maf/GDCdata", files.per.chunk = 100)
maf <- GDCprepare(q, directory = "data/tcga_maf/GDCdata")
saveRDS(as.data.frame(maf), "data/tcga_maf/tcga_brca_maf.rds")
cat("Variants:", nrow(maf), " Samples:", length(unique(maf$Tumor_Sample_Barcode)), "\nDone.\n")
