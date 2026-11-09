for (pk in c("survival", "glmnet")) {
  if (!requireNamespace(pk, quietly = TRUE)) stop("Package missing: ", pk)
}
suppressPackageStartupMessages({ library(survival); library(glmnet) })
set.seed(2026)
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

clin <- readRDS("data/clin_tumor_survival.rds")
expr <- readRDS("data/expr_tumor_vst.rds")
stopifnot(identical(rownames(clin), colnames(expr)))

## 1. Candidate genes = DE in BOTH paired and unpaired analyses
read_sig <- function(f) {
  d <- read.csv(f)
  d$gene_id[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) >= 1]
}
genes <- intersect(read_sig("results/DE_paired_tumor_vs_normal.csv"),
                   read_sig("results/DE_unpaired_tumor_vs_normal.csv"))
genes <- intersect(genes, rownames(expr))
cat("Overlap DE genes available in expression matrix:", length(genes), "\n")

## 2. 70:30 train/test split, stratified by event, fixed seed
train <- logical(nrow(clin))
for (e in c(0, 1)) {
  i <- which(clin$OS_event == e)
  train[sample(i, round(0.7 * length(i)))] <- TRUE
}
cat("Train:", sum(train), "patients,", sum(clin$OS_event[train]), "deaths\n")
cat("Test: ", sum(!train), "patients,", sum(clin$OS_event[!train]), "deaths\n")

## 3. Univariate Cox on TRAIN only
y_tr <- Surv(clin$OS_months[train], clin$OS_event[train])
X_tr <- t(expr[genes, train])
res <- t(sapply(seq_along(genes), function(j) {
  f <- tryCatch(coxph(y_tr ~ X_tr[, j]), error = function(e) NULL)
  if (is.null(f)) return(c(NA, NA))
  s <- summary(f)$coefficients
  c(s[1, "exp(coef)"], s[1, "Pr(>|z|)"])
}))
uni <- data.frame(gene_id = genes, HR = res[, 1], p = res[, 2])

ginfo <- tryCatch(as.data.frame(readRDS("data/gene_info.rds")), error = function(e) NULL)
add_name <- function(df) {
  if (!is.null(ginfo)) {
    id_col <- intersect(c("gene_id", "ensembl_gene_id", "id"), colnames(ginfo))[1]
    nm_col <- intersect(c("gene_name", "external_gene_name", "symbol"), colnames(ginfo))[1]
    if (!is.na(id_col) && !is.na(nm_col)) df$gene_name <- ginfo[[nm_col]][match(df$gene_id, ginfo[[id_col]])]
    else if (!is.na(nm_col)) df$gene_name <- ginfo[[nm_col]][match(df$gene_id, rownames(ginfo))]
  }
  df
}
write.csv(add_name(uni[order(uni$p), ]), "results/univariate_cox_train.csv", row.names = FALSE)
cand <- uni$gene_id[!is.na(uni$p) & uni$p < 0.05]
cat("Univariate Cox p<0.05 (train):", length(cand), "genes\n")

## 4. LASSO-Cox on TRAIN, max 10 genes (events-per-variable limit)
cv <- cv.glmnet(X_tr[, cand], y_tr, family = "cox", alpha = 1, nfolds = 10)
png("figures/lasso_cv.png", width = 1600, height = 1200, res = 220); plot(cv); dev.off()
ok <- which(cv$nzero >= 2 & cv$nzero <= 10)
if (length(ok) == 0) stop("No lambda gives 2-10 genes")
lam <- cv$lambda[ok[which.min(cv$cvm[ok])]]
co <- as.matrix(coef(cv$glmnet.fit, s = lam))
co <- co[co[, 1] != 0, , drop = FALSE]
sig <- add_name(data.frame(gene_id = rownames(co), coef = co[, 1]))
write.csv(sig, "results/signature_genes.csv", row.names = FALSE)
cat("\nSignature genes:", nrow(sig), "\n"); print(sig)

## 5. Risk score; cutoff = median of TRAIN (frozen)
rs <- as.numeric(t(expr[sig$gene_id, , drop = FALSE]) %*% sig$coef)
names(rs) <- colnames(expr)
cutoff <- median(rs[train])

report <- function(idx, label) {
  d <- data.frame(t = clin$OS_months[idx], e = clin$OS_event[idx], rs = rs[idx],
                  high = as.numeric(rs[idx] > cutoff))
  fit <- coxph(Surv(t, e) ~ high, data = d)
  ci  <- summary(fit)$conf.int
  lr  <- survdiff(Surv(t, e) ~ high, data = d)
  cidx <- concordance(Surv(t, e) ~ rs, data = d, reverse = TRUE)$concordance
  cat(sprintf("%s: n=%d deaths=%d | high-risk HR=%.2f (95%% CI %.2f-%.2f) | log-rank p=%.2e | C-index=%.3f\n",
              label, nrow(d), sum(d$e), ci[1, 1], ci[1, 3], ci[1, 4],
              1 - pchisq(lr$chisq, 1), cidx))
}
cat("\n")
report(which(train), "TRAIN")
report(which(!train), "TEST ")

saveRDS(list(genes = sig, cutoff = cutoff, lambda = lam, train_ids = colnames(expr)[train],
             test_ids = colnames(expr)[!train], risk = rs),
        "data/signature_model.rds")
cat("\nSaved: results/signature_genes.csv, data/signature_model.rds, figures/lasso_cv.png\n")
