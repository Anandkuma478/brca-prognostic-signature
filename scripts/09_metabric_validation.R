# 09_metabric_validation.R
# Independent validation of the frozen SCAN-B models (A and B) in METABRIC. No model is refitted.
# Primary endpoint: overall survival. Secondary endpoint: disease-specific survival.
# Run from the project root, once, after 06c_signature_check.R has finished:
#   Rscript scripts/09_metabric_validation.R 2>&1 | tee metabric.log
# Code test only (survival labels shuffled, so no real result is produced):
#   Rscript scripts/09_metabric_validation.R --shuffle <output_dir>

suppressPackageStartupMessages({ library(survival); library(data.table) })
set.seed(2030)

args    <- commandArgs(trailingOnly = TRUE)
shuffle <- "--shuffle" %in% args
extra   <- setdiff(args, "--shuffle")
outdir  <- if (length(extra)) extra[1] else "results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
models_file <- "data/signature_scanb_two_models.rds"
if (!shuffle && file.mtime(models_file) <= file.mtime("results/DE_paired_tumor_vs_normal.csv"))
  stop("Frozen models are older than the corrected DE results. Wait for 06c_signature_check.R to finish.")
if (shuffle) cat("*** SHUFFLED SURVIVAL LABELS: CODE TEST ONLY, NOT A RESULT ***\n")

read_cbio <- function(f) {
  L <- readLines(f, warn = FALSE)
  read.delim(text = L[!startsWith(L, "#")], quote = "", na.strings = c("", "NA"), check.names = FALSE)
}
pat <- read_cbio("data/metabric/raw/data_clinical_patient.txt")
sam <- read_cbio("data/metabric/raw/data_clinical_sample.txt")
models <- readRDS(models_file)
sigs   <- list(A = models$sigA, B = models$sigB)
hubs   <- read.csv("results/string_hub_genes.csv")
needed <- unique(c(unlist(lapply(sigs, function(s) s$symbol)), hubs$symbol[hubs$hub]))

## Expression: log2 intensity; duplicated symbols keep the highest mean
ex  <- fread("data/metabric/raw/data_mrna_illumina_microarray.txt", data.table = FALSE)
cat("METABRIC expression:", nrow(ex), "rows,", ncol(ex) - 2, "samples\n")
ex  <- ex[ex$Hugo_Symbol %in% needed, ]
mat <- as.matrix(ex[, -(1:2)]); rownames(mat) <- ex$Hugo_Symbol
if (anyDuplicated(rownames(mat))) {
  mat <- mat[order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
}

## Clinical data
d <- merge(pat, sam, by = "PATIENT_ID")
d <- d[d$PATIENT_ID %in% colnames(mat), ]
d$time <- as.numeric(d$OS_MONTHS)
d$os   <- ifelse(d$OS_STATUS == "1:DECEASED", 1, ifelse(d$OS_STATUS == "0:LIVING", 0, NA))
d$dss  <- ifelse(is.na(d$VITAL_STATUS), NA, as.numeric(d$VITAL_STATUS == "Died of Disease"))
n0 <- nrow(d)
d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$os), ]
cat("Patients with expression:", n0, " Excluded (missing or non-positive survival):", n0 - nrow(d),
    " Analysed:", nrow(d), "\n")
if (shuffle) { i <- sample(nrow(d)); d[, c("time", "os", "dss")] <- d[i, c("time", "os", "dss")] }

d$age    <- as.numeric(d$AGE_AT_DIAGNOSIS)
d$node   <- as.numeric(as.numeric(d$LYMPH_NODES_EXAMINED_POSITIVE) > 0)
d$grade  <- factor(d$GRADE)
d$er     <- ifelse(d$ER_IHC == "Positve", 1, ifelse(d$ER_IHC == "Negative", 0, NA))   # spelling as in source file
d$pam50  <- factor(ifelse(d$CLAUDIN_SUBTYPE == "NC", NA, d$CLAUDIN_SUBTYPE))
d$chemo  <- as.numeric(d$CHEMOTHERAPY == "YES")
d$endo   <- as.numeric(d$HORMONE_THERAPY == "YES")
if (!shuffle) {
  cat("Deaths:", sum(d$os), " Disease-specific deaths:", sum(d$dss, na.rm = TRUE), "\n")
  rf <- survfit(Surv(time, 1 - os) ~ 1, data = d)
  cat("Median follow-up (reverse KM, months):", round(summary(rf)$table["median"], 1), "\n")
  cat("Median age:", median(d$age, na.rm = TRUE), "\nPAM50:\n"); print(table(d$pam50, useNA = "ifany"))
  cat("ER (1=pos):\n"); print(table(d$er, useNA = "ifany")); cat("Node positive:\n"); print(table(d$node, useNA = "ifany"))
  cat("Grade:\n"); print(table(d$grade, useNA = "ifany")); cat("Chemotherapy:\n"); print(table(d$chemo, useNA = "ifany"))
  cat("Endocrine therapy:\n"); print(table(d$endo, useNA = "ifany"))
}

## Standardise within METABRIC; a gene absent or unmeasured contributes zero
z <- t(scale(t(mat[, d$PATIENT_ID, drop = FALSE])))
z[is.na(z)] <- 0

cidx <- function(y, lp) concordance(y ~ lp, reverse = TRUE)$concordance
clin_terms <- "age + node + grade + er + pam50 + chemo + endo"

evaluate <- function(nm, endpoint) {
  s  <- sigs[[nm]]
  in_data <- s$symbol %in% rownames(z)
  rs <- as.numeric(t(z[s$symbol[in_data], , drop = FALSE]) %*% s$coef[in_data])
  keep <- !is.na(d[[endpoint]])
  x <- data.frame(t = d$time, e = d[[endpoint]], score = as.numeric(scale(rs)),
                  d[, c("age", "node", "grade", "er", "pam50", "chemo", "endo")])[keep, ]
  y  <- Surv(x$t, x$e)
  u  <- summary(coxph(y ~ x$score))$conf.int[1, ]
  dd <- na.omit(x)
  m0 <- coxph(as.formula(paste("Surv(t, e) ~", clin_terms)), data = dd)
  m1 <- coxph(as.formula(paste("Surv(t, e) ~ score +", clin_terms)), data = dd)
  a  <- summary(m1)$conf.int["score", ]
  data.frame(model = nm, endpoint = endpoint, genes_in_model = nrow(s), genes_matched = sum(in_data),
             n = nrow(x), events = sum(x$e), c_index = cidx(y, x$score),
             HR_per_SD = u[1], HR_low = u[3], HR_high = u[4],
             n_adjusted = nrow(dd), events_adjusted = sum(dd$e),
             adj_HR = a[1], adj_low = a[3], adj_high = a[4],
             adj_p = summary(m1)$coefficients["score", "Pr(>|z|)"],
             LR_p_vs_clinical = anova(m0, m1)$`Pr(>|Chi|)`[2],
             c_clinical = concordance(m0)$concordance, c_clinical_score = concordance(m1)$concordance,
             row.names = NULL)
}
out <- do.call(rbind, lapply(c("os", "dss"), function(e) do.call(rbind, lapply(c("A", "B"), evaluate, endpoint = e))))
print(format(out, digits = 3), row.names = FALSE)
write.csv(out, file.path(outdir, if (shuffle) "metabric_SHUFFLED_TEST.csv" else "metabric_validation.csv"), row.names = FALSE)
if (!shuffle) saveRDS(list(clin = d, z = z), "data/metabric_prepared.rds")  # read by 12_survival_evaluation.R
cat("\nDone.\n")
