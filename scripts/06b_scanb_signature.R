suppressPackageStartupMessages({ library(survival); library(glmnet); library(doParallel) }); registerDoParallel(cores = 8)
set.seed(2026)
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

sc <- readRDS("data/scanb_clin.rds");        sx <- readRDS("data/scanb_expr.rds")
tc <- readRDS("data/clin_tumor_survival.rds"); tx <- readRDS("data/expr_tumor_vst.rds")
stopifnot(identical(rownames(sc), colnames(sx)), identical(rownames(tc), colnames(tx)))

## gene map (TCGA Ensembl id -> symbol, type)
ginfo <- as.data.frame(readRDS("data/gene_info.rds"))
id_col   <- intersect(c("gene_id", "ensembl_gene_id", "id"), colnames(ginfo))[1]
nm_col   <- intersect(c("gene_name", "external_gene_name", "symbol"), colnames(ginfo))[1]
type_col <- intersect(c("gene_type", "gene_biotype", "biotype"), colnames(ginfo))[1]
stopifnot(!is.na(nm_col))
map <- data.frame(gene_id = if (!is.na(id_col)) ginfo[[id_col]] else rownames(ginfo),
                  symbol = ginfo[[nm_col]],
                  type = if (!is.na(type_col)) ginfo[[type_col]] else NA, stringsAsFactors = FALSE)

## 1. Candidates (all filters are outcome-blind)
read_sig <- function(f) { d <- read.csv(f); d$gene_id[!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) >= 1] }
de <- intersect(read_sig("results/DE_paired_tumor_vs_normal.csv"), read_sig("results/DE_unpaired_tumor_vs_normal.csv"))
cand <- map[map$gene_id %in% de & map$gene_id %in% rownames(tx), ]
cat("DE in both analyses (in TCGA matrix):", nrow(cand), "\n")
if (!all(is.na(cand$type))) cand <- cand[cand$type == "protein_coding", ]
cat("Protein-coding:", nrow(cand), "\n")
cand <- cand[!is.na(cand$symbol) & !duplicated(cand$symbol) & cand$symbol %in% rownames(sx), ]
cat("Present in SCAN-B (unique symbols):", nrow(cand), "\n")
det  <- rowMeans(sx[cand$symbol, ] > -3.3) >= 0.5
tok  <- rowMeans(tx[cand$gene_id, ]) >= median(rowMeans(tx))
cand <- cand[det & tok, ]
cat("Expressed in both cohorts:", nrow(cand), "\n")
v <- apply(sx[cand$symbol, ], 1, var)
cand <- cand[order(v, decreasing = TRUE), ][seq_len(min(2000, nrow(cand))), ]
cat("Candidates for LASSO:", nrow(cand), "\n")

Xs <- scale(t(sx[cand$symbol, ]))                       # SCAN-B, z per gene
Xt <- scale(t(tx[cand$gene_id, ])); colnames(Xt) <- cand$symbol   # TCGA, z per gene
ys <- Surv(sc$OS_months, sc$OS_event); yt <- Surv(tc$OS_months, tc$OS_event)
cidx <- function(y, lp) if (sd(lp) == 0) 0.5 else concordance(y ~ lp, reverse = TRUE)$concordance

## 2. Nested CV in SCAN-B
K <- 10; R <- 3; n <- nrow(Xs)
lp1 <- lpm <- matrix(0, n, R)
for (r in seq_len(R)) {
  folds <- integer(n)
  for (e in c(0, 1)) { i <- which(sc$OS_event == e); folds[i] <- sample(rep(seq_len(K), length.out = length(i))) }
  for (k in seq_len(K)) {
    tr <- folds != k
    cvk <- cv.glmnet(Xs[tr, ], ys[tr], family = "cox", alpha = 1, nfolds = 5, parallel = TRUE)
    lp1[!tr, r] <- as.numeric(predict(cvk, Xs[!tr, ], s = "lambda.1se", type = "link"))
    lpm[!tr, r] <- as.numeric(predict(cvk, Xs[!tr, ], s = "lambda.min", type = "link"))
  }
  cat("CV repeat", r, "done\n")
}
c1 <- sapply(seq_len(R), function(r) cidx(ys, lp1[, r])); cm <- sapply(seq_len(R), function(r) cidx(ys, lpm[, r]))
cat("\nSCAN-B nested-CV C-index, lambda.1se:", round(c1, 3), "| mean", round(mean(c1), 3), "\n")
cat("SCAN-B nested-CV C-index, lambda.min:", round(cm, 3), "| mean", round(mean(cm), 3), "\n")

## 3. Final model on all SCAN-B
cvf <- cv.glmnet(Xs, ys, family = "cox", alpha = 1, nfolds = 10, parallel = TRUE)
png("figures/lasso_cv_scanb.png", width = 1600, height = 1200, res = 220); plot(cvf); dev.off()
n1 <- sum(as.matrix(coef(cvf, s = "lambda.1se")) != 0); nm <- sum(as.matrix(coef(cvf, s = "lambda.min")) != 0)
use <- if (n1 >= 3) "lambda.1se" else "lambda.min"
cat("Nonzero genes: lambda.1se =", n1, " lambda.min =", nm, " -> using", use, "\n")
co <- as.matrix(coef(cvf, s = use)); co <- co[co[, 1] != 0, , drop = FALSE]
sig <- data.frame(symbol = rownames(co), coef = co[, 1])
sig$gene_id <- cand$gene_id[match(sig$symbol, cand$symbol)]
write.csv(sig, "results/signature_scanb.csv", row.names = FALSE)
cat("\nSignature genes:", nrow(sig), "\n"); print(sig)

rs_s <- as.numeric(Xs[, sig$symbol, drop = FALSE] %*% sig$coef)
rs_t <- as.numeric(Xt[, sig$symbol, drop = FALSE] %*% sig$coef)
cut_s <- median(rs_s)                                    # frozen on SCAN-B

km <- function(y, d, grp, file, ttl) {
  png(file, width = 1600, height = 1200, res = 220)
  f <- survfit(y ~ grp); plot(f, col = c("blue", "red"), lwd = 2, xlab = "Months", ylab = "Overall survival", main = ttl)
  legend("bottomleft", c("Low risk", "High risk"), col = c("blue", "red"), lwd = 2, bty = "n"); dev.off()
}
evalc <- function(y, rs, label, file) {
  hi <- as.numeric(rs > cut_s)
  z <- as.numeric(scale(rs))
  h1 <- summary(coxph(y ~ z))$conf.int[1, ]; h2 <- summary(coxph(y ~ hi))$conf.int[1, ]
  lr <- survdiff(y ~ hi)
  cat(sprintf("%s: C-index=%.3f | per-SD HR=%.2f (%.2f-%.2f) | high-vs-low (SCAN-B cutoff) HR=%.2f (%.2f-%.2f), log-rank p=%.3g, n high=%d\n",
              label, cidx(y, rs), h1[1], h1[3], h1[4], h2[1], h2[3], h2[4], 1 - pchisq(lr$chisq, 1), sum(hi)))
  if (length(unique(hi)) == 2) km(y, NULL, hi, file, label)
}
cat("\n")
evalc(ys, rs_s, "SCAN-B (apparent, optimistic)", "figures/km_scanb_apparent.png")
evalc(yt, rs_t, "TCGA (validation)", "figures/km_tcga_validation.png")

## 4. Independence from clinical variables
d1 <- data.frame(t = sc$OS_months, e = sc$OS_event, lp = as.numeric(scale(rowMeans(lp1))),
                 age = sc$age, node = sc$node_pos, grade = sc$nhg, er = sc$er, pam50 = sc$pam50,
                 chemo = sc$chemo, endo = sc$endocrine)
m1 <- coxph(Surv(t, e) ~ lp + age + node + grade + er + pam50 + chemo + endo, data = d1)
s1 <- summary(m1)$conf.int["lp", ]
cat(sprintf("\nSCAN-B multivariate (CV risk score per SD): HR=%.2f (%.2f-%.2f), p=%.3g, n=%d\n",
            s1[1], s1[3], s1[4], summary(m1)$coefficients["lp", "Pr(>|z|)"], m1$n))
d2 <- data.frame(t = tc$OS_months, e = tc$OS_event, lp = as.numeric(scale(rs_t)), age = tc$age, stage = tc$stage, pam50 = tc$pam50)
m2 <- coxph(Surv(t, e) ~ lp + age + stage + pam50, data = d2)
s2 <- summary(m2)$conf.int["lp", ]
cat(sprintf("TCGA multivariate (risk score per SD): HR=%.2f (%.2f-%.2f), p=%.3g, n=%d\n",
            s2[1], s2[3], s2[4], summary(m2)$coefficients["lp", "Pr(>|z|)"], m2$n))

saveRDS(list(sig = sig, cutoff = cut_s, used = use, candidates = cand, cv_c1 = c1, cv_cm = cm),
        "data/signature_scanb_model.rds")
cat("\nSaved: results/signature_scanb.csv, data/signature_scanb_model.rds, figures/*.png\n")
