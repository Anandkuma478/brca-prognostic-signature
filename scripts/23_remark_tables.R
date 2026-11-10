# 23_remark_tables.R
# Two tables requested by the REMARK reporting recommendations, from the frozen scores (nothing is refitted):
#   1. Full multivariable Cox models: hazard ratios for the risk score and for every clinical variable
#      (same complete-case models as 06c_signature_check.R and 09_metabric_validation.R).
#   2. Relation of the risk score to standard prognostic variables.
# Run from the project root with the brca_genomics environment, after 14_cindex_ci.R:
#   Rscript scripts/23_remark_tables.R 2>&1 | tee remark_tables.log

suppressPackageStartupMessages(library(survival))

models <- readRDS("data/signature_scanb_two_models.rds"); sigs <- list(A = models$sigA, B = models$sigB)
cv <- readRDS("data/nested_cv_predictions.rds")
gi <- as.data.frame(readRDS("data/gene_info.rds")); sym2id <- setNames(gi$gene_id, gi$gene_name)
zscore <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
score <- function(z, s) { k <- s$symbol %in% rownames(z); as.numeric(t(z[s$symbol[k], , drop = FALSE]) %*% s$coef[k]) }

sc <- readRDS("data/scanb_clin.rds")
tc <- readRDS("data/clin_tumor_survival.rds"); tx <- readRDS("data/expr_tumor_vst.rds")
tg <- unique(c(sigs$A$symbol, sigs$B$symbol)); id <- sym2id[tg]; keep <- !is.na(id) & id %in% rownames(tx)
zt <- zscore(tx[id[keep], ]); rownames(zt) <- tg[keep]
mb <- readRDS("data/metabric_prepared.rds"); mc <- mb$clin

cohort_data <- function(m) list(
  `SCAN-B` = data.frame(t = sc$OS_months, e = sc$OS_event, score = as.numeric(scale(rowMeans(if (m == "A") cv$lpA else cv$lpB))),
                        age = sc$age, node = sc$node_pos, grade = sc$nhg, er = sc$er, pam50 = sc$pam50, chemo = sc$chemo, endo = sc$endocrine),
  `TCGA-BRCA` = data.frame(t = tc$OS_months, e = tc$OS_event, score = as.numeric(scale(score(zt, sigs[[m]]))),
                           age = tc$age, stage = tc$stage, pam50 = tc$pam50),
  METABRIC = data.frame(t = mc$time, e = mc$os, score = as.numeric(scale(score(mb$z, sigs[[m]]))),
                        mc[, c("age", "node", "grade", "er", "pam50", "chemo", "endo")]))

## 1. Full multivariable models
full <- do.call(rbind, lapply(c("A", "B"), function(m) {
  dl <- cohort_data(m)
  do.call(rbind, lapply(names(dl), function(co) {
    dd <- na.omit(dl[[co]])
    fit <- coxph(Surv(t, e) ~ ., data = dd)
    s <- summary(fit)
    data.frame(cohort = co, model = m, n = nrow(dd), events = sum(dd$e), term = rownames(s$conf.int),
               HR = s$conf.int[, 1], HR_low = s$conf.int[, 3], HR_high = s$conf.int[, 4],
               p = s$coefficients[, "Pr(>|z|)"], row.names = NULL)
  }))
}))
## The score rows must reproduce the adjusted hazard ratios reported in the manuscript
chk <- full[full$term == "score", ]
expect <- c(`SCAN-B A` = 1.38, `SCAN-B B` = 1.47, `TCGA-BRCA A` = 1.50, `TCGA-BRCA B` = 1.16, `METABRIC A` = 1.17, `METABRIC B` = 1.12)
got <- setNames(round(chk$HR, 2), paste(chk$cohort, chk$model))
print(rbind(expected = expect[names(got)], reproduced = got))
stopifnot(isTRUE(all.equal(unname(expect[names(got)]), unname(got))))
write.csv(full, "results/multivariable_models_full.csv", row.names = FALSE)

## 2. Risk score (Model A, standardised) by standard prognostic variable
rel <- do.call(rbind, lapply(names(cohort_data("A")), function(co) {
  d <- cohort_data("A")[[co]]
  vars <- intersect(c("node", "grade", "stage", "er", "pam50"), names(d))
  out <- lapply(vars, function(v) {
    x <- d[!is.na(d[[v]]), ]; g <- factor(x[[v]])
    p <- kruskal.test(x$score ~ g)$p.value
    do.call(rbind, lapply(levels(g), function(l) {
      s <- x$score[g == l]
      data.frame(cohort = co, variable = v, level = l, n = length(s), median = median(s),
                 q1 = unname(quantile(s, 0.25)), q3 = unname(quantile(s, 0.75)), p_kruskal = p)
    }))
  })
  a <- d[!is.na(d$age), ]; ct <- cor.test(a$age, a$score, method = "spearman", exact = FALSE)
  rbind(do.call(rbind, out), data.frame(cohort = co, variable = "age", level = "Spearman rho", n = nrow(a),
                                        median = unname(ct$estimate), q1 = NA, q3 = NA, p_kruskal = ct$p.value))
}))
write.csv(rel, "results/score_vs_clinical.csv", row.names = FALSE)
print(format(rel, digits = 3), row.names = FALSE)
cat("\nSaved: results/multivariable_models_full.csv, results/score_vs_clinical.csv\n")
