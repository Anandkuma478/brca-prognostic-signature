# 14_cindex_ci.R
# 95% confidence intervals for the C-index of the frozen models, and a proportional-hazards check of each score.
# Uses the scores of the frozen models only; nothing is refitted. Run after 09_metabric_validation.R:
#   Rscript scripts/14_cindex_ci.R 2>&1 | tee cindex_ci.log

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
mb <- readRDS("data/metabric_prepared.rds")

one <- function(cohort, endpoint, m, time, event, lp) {
  ok <- !is.na(event)
  y  <- Surv(time[ok], event[ok]); lp <- lp[ok]
  co <- concordance(y ~ lp, reverse = TRUE); se <- sqrt(co$var)
  ph <- cox.zph(coxph(y ~ lp))$table["lp", "p"]
  sdz <- as.numeric(scale(lp)); h <- summary(coxph(y ~ sdz))$conf.int[1, ]
  data.frame(cohort = cohort, endpoint = endpoint, model = m, n = length(lp), events = sum(event[ok]),
             c_index = co$concordance, c_low = co$concordance - 1.96 * se, c_high = co$concordance + 1.96 * se,
             HR_per_SD = h[1], HR_low = h[3], HR_high = h[4], PH_test_p = ph, row.names = NULL)
}
out <- do.call(rbind, lapply(c("A", "B"), function(m) rbind(
  one("SCAN-B (mean out-of-fold)", "OS", m, sc$OS_months, sc$OS_event, rowMeans(if (m == "A") cv$lpA else cv$lpB)),
  one("TCGA-BRCA", "OS", m, tc$OS_months, tc$OS_event, score(zt, sigs[[m]])),
  one("METABRIC", "OS", m, mb$clin$time, mb$clin$os, score(mb$z, sigs[[m]])),
  one("METABRIC", "DSS", m, mb$clin$time, mb$clin$dss, score(mb$z, sigs[[m]])))))
print(format(out, digits = 3), row.names = FALSE)
write.csv(out, "results/cindex_summary.csv", row.names = FALSE)

## TCGA-BRCA: added value over age, stage and PAM50 (complete cases)
added <- do.call(rbind, lapply(c("A", "B"), function(m) {
  d <- na.omit(data.frame(t = tc$OS_months, e = tc$OS_event, lp = as.numeric(scale(score(zt, sigs[[m]]))), age = tc$age, stage = tc$stage, pam50 = tc$pam50))
  m0 <- coxph(Surv(t, e) ~ age + stage + pam50, data = d); m1 <- coxph(Surv(t, e) ~ lp + age + stage + pam50, data = d)
  data.frame(model = m, n = nrow(d), events = sum(d$e), LR_p_vs_clinical = anova(m0, m1)$`Pr(>|Chi|)`[2],
             c_clinical = concordance(m0)$concordance, c_clinical_score = concordance(m1)$concordance)
}))
print(format(added, digits = 3), row.names = FALSE)
write.csv(added, "results/tcga_added_value.csv", row.names = FALSE)
cat("\nSaved: results/cindex_summary.csv, results/tcga_added_value.csv\n")
