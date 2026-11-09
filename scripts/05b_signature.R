suppressPackageStartupMessages({ library(survival); library(glmnet) })
set.seed(2026)
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

clin <- readRDS("data/clin_tumor_survival.rds")
expr <- readRDS("data/expr_tumor_vst.rds")
stopifnot(identical(rownames(clin), colnames(expr)))
y <- Surv(clin$OS_months, clin$OS_event)

## 1. Candidates: DE in both analyses, then outcome-blind filters
read_sig <- function(f) {
  d <- read.csv(f)
  d$gene_id[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) >= 1]
}
genes <- intersect(read_sig("results/DE_paired_tumor_vs_normal.csv"),
                   read_sig("results/DE_unpaired_tumor_vs_normal.csv"))
genes <- intersect(genes, rownames(expr))
cat("DE in both analyses:", length(genes), "\n")

ginfo <- tryCatch(as.data.frame(readRDS("data/gene_info.rds")), error = function(e) NULL)
id_col <- nm_col <- type_col <- NA
if (!is.null(ginfo)) {
  id_col   <- intersect(c("gene_id", "ensembl_gene_id", "id"), colnames(ginfo))[1]
  nm_col   <- intersect(c("gene_name", "external_gene_name", "symbol"), colnames(ginfo))[1]
  type_col <- intersect(c("gene_type", "gene_biotype", "biotype"), colnames(ginfo))[1]
}
ginfo_ids <- if (!is.null(ginfo)) { if (!is.na(id_col)) ginfo[[id_col]] else rownames(ginfo) } else NULL
if (!is.null(ginfo) && !is.na(type_col)) {
  pc <- ginfo_ids[ginfo[[type_col]] == "protein_coding"]
  genes <- genes[genes %in% pc]
  cat("After protein-coding filter:", length(genes), "\n")
} else {
  cat("No gene type column in gene_info.rds: protein-coding filter SKIPPED\n")
}
floor_mu <- median(rowMeans(expr))
genes <- genes[rowMeans(expr[genes, ]) >= floor_mu]
cat("After expression floor:", length(genes), "\n")
v <- apply(expr[genes, ], 1, var)
genes <- names(sort(v, decreasing = TRUE))[seq_len(min(2000, length(v)))]
cat("Candidates used in LASSO (top variance):", length(genes), "\n")
X <- t(expr[genes, ])

add_name <- function(df) {
  if (!is.null(ginfo) && !is.na(nm_col)) df$gene_name <- ginfo[[nm_col]][match(df$gene_id, ginfo_ids)]
  df
}

## 2. Nested cross-validation (selection + LASSO inside every training fold)
K <- 10; R <- 3
lp_cv <- matrix(0, nrow(X), R)
for (r in seq_len(R)) {
  folds <- integer(nrow(X))
  for (e in c(0, 1)) {
    i <- which(clin$OS_event == e)
    folds[i] <- sample(rep(seq_len(K), length.out = length(i)))
  }
  for (k in seq_len(K)) {
    tr <- folds != k
    cvk <- cv.glmnet(X[tr, ], y[tr], family = "cox", alpha = 1, nfolds = 5)
    lp_cv[!tr, r] <- as.numeric(predict(cvk, X[!tr, ], s = "lambda.1se", type = "link"))
  }
  cat("CV repeat", r, "done\n")
}
cidx <- sapply(seq_len(R), function(r) concordance(y ~ lp_cv[, r], reverse = TRUE)$concordance)
cat("\nCross-validated C-index per repeat:", round(cidx, 3), "| mean =", round(mean(cidx), 3), "\n")

lp <- rowMeans(lp_cv)
if (sd(lp) > 0) {
  high <- as.numeric(lp > median(lp))
  fit <- coxph(y ~ high); ci <- summary(fit)$conf.int
  lr <- survdiff(y ~ high)
  cat(sprintf("CV risk groups: HR=%.2f (95%% CI %.2f-%.2f), log-rank p=%.3g\n",
              ci[1, 1], ci[1, 3], ci[1, 4], 1 - pchisq(lr$chisq, 1)))
  d <- data.frame(OS_months = clin$OS_months, OS_event = clin$OS_event,
                  lp = as.numeric(scale(lp)), age = clin$age, stage = clin$stage, pam50 = clin$pam50)
  for (f in c("lp + age + pam50", "lp + age + stage + pam50")) {
    m <- coxph(as.formula(paste("Surv(OS_months, OS_event) ~", f)), data = d)
    s <- summary(m)$conf.int
    cat(sprintf("Multivariate (%s), n=%d: CV risk score per SD HR=%.2f (95%% CI %.2f-%.2f), p=%.3g\n",
                f, m$n, s["lp", 1], s["lp", 3], s["lp", 4], summary(m)$coefficients["lp", "Pr(>|z|)"]))
  }
} else {
  cat("CV predictions are constant: LASSO selected no genes in the folds (no usable signal).\n")
}

## 3. Final model on all patients
cvf <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)
png("figures/lasso_cv_final.png", width = 1600, height = 1200, res = 220); plot(cvf); dev.off()
for (s in c("lambda.1se", "lambda.min")) {
  co <- as.matrix(coef(cvf, s = s)); cat(s, ": nonzero genes =", sum(co[, 1] != 0), "\n")
}
co <- as.matrix(coef(cvf, s = "lambda.1se")); co <- co[co[, 1] != 0, , drop = FALSE]
if (nrow(co) > 0) {
  sig <- add_name(data.frame(gene_id = rownames(co), coef = co[, 1]))
  write.csv(sig, "results/signature_genes_final.csv", row.names = FALSE)
  print(sig)
  rs <- as.numeric(X[, sig$gene_id, drop = FALSE] %*% sig$coef)
  cat("Apparent (optimistic) C-index of final model:",
      round(concordance(y ~ rs, reverse = TRUE)$concordance, 3), "\n")
  saveRDS(list(genes = sig, cutoff = median(rs), lambda = cvf$lambda.1se, candidates = genes),
          "data/signature_model_final.rds")
  cat("Saved: results/signature_genes_final.csv, data/signature_model_final.rds\n")
} else {
  cat("\nNo genes selected at lambda.1se: no stable signature in TCGA-BRCA with these data.\n")
}
