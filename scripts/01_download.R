library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("/home/administrator/Desktop/anand 2/brca genomics")

query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, files.per.chunk = 20)
se <- GDCprepare(query)
saveRDS(se, "data/tcga_brca_se.rds")

clin <- GDCquery_clinic("TCGA-BRCA", type = "clinical")
saveRDS(clin, "data/tcga_brca_clinical.rds")

cat("Samples:", ncol(se), " Genes:", nrow(se), "\n")
