suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})

counts <- readRDS("data/counts_filtered.rds")
if (inherits(counts, "SummarizedExperiment")) counts <- assay(counts, "unstranded")
counts <- as.matrix(counts)
meta <- as.data.frame(readRDS("data/meta_clean.rds"))
meta <- meta[colnames(counts), , drop = FALSE]

bc <- colnames(counts)
is_tumor <- substr(bc, 14, 15) == "01"
meta <- meta[is_tumor, ]
counts <- counts[, is_tumor]
meta$patient <- substr(rownames(meta), 1, 12)
stopifnot(!any(duplicated(meta$patient)))
cat("Primary tumors:", nrow(meta), "\n")

dtd <- as.numeric(meta$days_to_death)
dlf <- as.numeric(meta$days_to_last_follow_up)
meta$OS_event <- ifelse(meta$vital_status == "Dead", 1, ifelse(meta$vital_status == "Alive", 0, NA))
meta$OS_time  <- ifelse(meta$OS_event == 1, dtd, dlf)

meta$age <- as.numeric(meta$age_at_index)
st <- gsub("Stage ", "", as.character(meta$ajcc_pathologic_stage))
st <- sub("[ABC]$", "", st)
st[!st %in% c("I", "II", "III", "IV")] <- NA
meta$stage <- factor(st, levels = c("I", "II", "III", "IV"))
meta$pam50 <- factor(meta$paper_BRCA_Subtype_PAM50)

bad <- is.na(meta$OS_time) | is.na(meta$OS_event) | meta$OS_time <= 0
cat("Removed (NA / time<=0):", sum(bad), "\n")
meta <- meta[!bad, ]
counts <- counts[, rownames(meta)]
meta$OS_months <- meta$OS_time / 30.44

cat("\nFinal tumors:", nrow(meta), "\n")
cat("Deaths:", sum(meta$OS_event), " Alive/censored:", sum(meta$OS_event == 0), "\n")
cat("Median follow-up (months):", round(median(meta$OS_months), 1), "\n")
cat("Stage:\n"); print(table(meta$stage, useNA = "ifany"))
cat("PAM50:\n"); print(table(meta$pam50, useNA = "ifany"))

dds <- DESeqDataSetFromMatrix(round(counts), colData = meta, design = ~ 1)
vsd <- vst(dds, blind = TRUE)
expr <- assay(vsd)

saveRDS(meta, "data/clin_tumor_survival.rds")
saveRDS(expr, "data/expr_tumor_vst.rds")
cat("\nSaved: data/clin_tumor_survival.rds, data/expr_tumor_vst.rds\n")
