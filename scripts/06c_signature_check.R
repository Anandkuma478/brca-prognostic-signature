suppressPackageStartupMessages({ library(survival); library(glmnet); library(doParallel) })
registerDoParallel(cores = 8)
set.seed(2027)
sc <- readRDS("data/scanb_clin.rds");          sx <- readRDS("data/scanb_expr.rds")
tc <- readRDS("data/clin_tumor_survival.rds"); tx <- readRDS("data/expr_tumor_vst.rds")
cand <- readRDS("data/signature_scanb_model.rds")$candidates
Xs <- scale(t(sx[cand$symbol, ]))
Xt <- scale(t(tx[cand$gene_id, ])); colnames(Xt) <- cand$symbol
ys <- Surv(sc$OS_months, sc$OS_event); yt <- Surv(tc$OS_months, tc$OS_event)
cidx <- function(y, lp) if (sd(lp) == 0) 0.5 else concordance(y ~ lp, reverse = TRUE)$concordance
sparse_lam <- function(cv, maxg = 30) min(cv$lambda[cv$nzero <= maxg])

## 1. Nested CV: A = lambda.min (full), B = sparse (<= 30 genes)
K <- 5; R <- 2; n <- nrow(Xs)
lpA <- lpB <- matrix(0, n, R)
for (r in seq_len(R)) {
  folds <- integer(n)
  for (e in c(0, 1)) { i <- which(sc$OS_event == e); folds[i] <- sample(rep(seq_len(K), length.out = length(i))) }
  for (k in seq_len(K)) {
    tr <- folds != k
    cvk <- cv.glmnet(Xs[tr, ], ys[tr], family = "cox", alpha = 1, nfolds = 5, parallel = TRUE)
    lpA[!tr, r] <- as.numeric(predict(cvk, Xs[!tr, ], s = "lambda.min", type = "link"))
    lpB[!tr, r] <- as.numeric(predict(cvk, Xs[!tr, ], s = sparse_lam(cvk), type = "link"))
    cat("repeat", r, "fold", k, "done\n")
  }
}
saveRDS(list(lpA = lpA, lpB = lpB), "data/nested_cv_predictions.rds")
cA <- sapply(seq_len(R), function(r) cidx(ys, lpA[, r])); cB <- sapply(seq_len(R), function(r) cidx(ys, lpB[, r]))
cat("\nSCAN-B nested-CV C-index  A (full, lambda.min):", round(cA, 3), "| mean", round(mean(cA), 3), "\n")
cat("SCAN-B nested-CV C-index  B (sparse <=30 genes):", round(cB, 3), "| mean", round(mean(cB), 3), "\n")
er <- sc$er
for (g in c(1, 0)) {
  idx <- which(er == g)
  cat(sprintf("ER=%d  n=%d deaths=%d | CV C-index A=%.3f  B=%.3f\n", g, length(idx), sum(sc$OS_event[idx]),
              cidx(ys[idx], rowMeans(lpA)[idx]), cidx(ys[idx], rowMeans(lpB)[idx])))
}

## 2. Final models on all SCAN-B
cvf <- cv.glmnet(Xs, ys, family = "cox", alpha = 1, nfolds = 10, parallel = TRUE)
getsig <- function(s) { co <- as.matrix(coef(cvf, s = s)); co <- co[co[, 1] != 0, , drop = FALSE]
                        data.frame(symbol = rownames(co), coef = co[, 1]) }
sigA <- getsig("lambda.min"); sigB <- getsig(sparse_lam(cvf))
write.csv(sigB, "results/signature_scanb_sparse30.csv", row.names = FALSE)
cat("\nFinal sizes: A =", nrow(sigA), "genes, B =", nrow(sigB), "genes\n"); print(sigB)

## 3. TCGA validation (frozen)
d2 <- data.frame(t = tc$OS_months, e = tc$OS_event, age = tc$age, stage = tc$stage, pam50 = tc$pam50)
for (nm in c("A", "B")) {
  sig <- if (nm == "A") sigA else sigB
  rs <- as.numeric(Xt[, sig$symbol, drop = FALSE] %*% sig$coef)
  z <- as.numeric(scale(rs)); h <- summary(coxph(yt ~ z))$conf.int[1, ]
  d2$lp <- z; m <- coxph(Surv(t, e) ~ lp + age + stage + pam50, data = d2); s <- summary(m)$conf.int["lp", ]
  cat(sprintf("TCGA model %s: C-index=%.3f | per-SD HR=%.2f (%.2f-%.2f) | adjusted HR=%.2f (%.2f-%.2f), p=%.3g, n=%d\n",
              nm, cidx(yt, rs), h[1], h[3], h[4], s[1], s[3], s[4], summary(m)$coefficients["lp", "Pr(>|z|)"], m$n))
}

## 4. SCAN-B independence (CORRECT CV scores) + added value over clinical model
d1 <- data.frame(t = sc$OS_months, e = sc$OS_event, age = sc$age, node = sc$node_pos, grade = sc$nhg,
                 er = sc$er, pam50 = sc$pam50, chemo = sc$chemo, endo = sc$endocrine)
for (nm in c("A", "B")) {
  d1$lp <- as.numeric(scale(rowMeans(if (nm == "A") lpA else lpB)))
  dd <- na.omit(d1)
  m0 <- coxph(Surv(t, e) ~ age + node + grade + er + pam50 + chemo + endo, data = dd)
  m1 <- coxph(Surv(t, e) ~ lp + age + node + grade + er + pam50 + chemo + endo, data = dd)
  s <- summary(m1)$conf.int["lp", ]
  cat(sprintf("SCAN-B model %s (n=%d): CV score per-SD adjusted HR=%.2f (%.2f-%.2f), p=%.3g | LR test vs clinical-only p=%.3g | C-index clinical=%.3f, clinical+score=%.3f\n",
              nm, nrow(dd), s[1], s[3], s[4], summary(m1)$coefficients["lp", "Pr(>|z|)"],
              anova(m0, m1)$`Pr(>|Chi|)`[2], concordance(m0)$concordance, concordance(m1)$concordance))
}
saveRDS(list(sigA = sigA, sigB = sigB), "data/signature_scanb_two_models.rds")
cat("\nSaved: results/signature_scanb_sparse30.csv, data/nested_cv_predictions.rds\n")
