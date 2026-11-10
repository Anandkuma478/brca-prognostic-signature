# 12_survival_evaluation.R
# Phase 3: Kaplan-Meier curves, time-dependent ROC (1, 3 and 5 years) and calibration for Models A and B,
# and univariate Cox analysis of the hub genes. Nothing is refitted on validation cohorts.
# SCAN-B uses out-of-fold scores from nested cross-validation; TCGA-BRCA and METABRIC use the frozen models.
# METABRIC is included only after 09_metabric_validation.R has produced data/metabric_prepared.rds.
# Run from the project root with the brca_modules environment, after 06c and 09:
#   Rscript scripts/12_survival_evaluation.R 2>&1 | tee survival_evaluation.log
# Code test only (survival labels shuffled, METABRIC never read):
#   Rscript scripts/12_survival_evaluation.R --shuffle <output_dir>

suppressPackageStartupMessages({ library(survival); library(timeROC); library(survminer); library(ggplot2) })
set.seed(2033)

args    <- commandArgs(trailingOnly = TRUE)
shuffle <- "--shuffle" %in% args
extra   <- setdiff(args, "--shuffle")
outdir  <- if (length(extra)) extra[1] else "."
resdir  <- file.path(outdir, "results"); figdir <- file.path(outdir, "figures")
dir.create(resdir, showWarnings = FALSE, recursive = TRUE); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
if (!shuffle && file.mtime("data/signature_scanb_two_models.rds") <= file.mtime("results/DE_paired_tumor_vs_normal.csv"))
  stop("Frozen models are older than the corrected DE results. Wait for 06c_signature_check.R to finish.")
if (shuffle) cat("*** SHUFFLED SURVIVAL LABELS: CODE TEST ONLY, NOT A RESULT ***\n")

TIMES <- c(12, 36, 60)   # months
CAL_T <- 36              # calibration horizon, months

models <- readRDS("data/signature_scanb_two_models.rds"); sigs <- list(A = models$sigA, B = models$sigB)
cand   <- readRDS("data/signature_scanb_model.rds")$candidates
cv     <- readRDS("data/nested_cv_predictions.rds")
hubs   <- read.csv("results/string_hub_genes.csv"); hubs <- hubs[hubs$hub, ]
gi     <- as.data.frame(readRDS("data/gene_info.rds"))

zscore <- function(m) { z <- t(scale(t(m))); z[is.na(z)] <- 0; z }
score  <- function(z, s) as.numeric(t(z[s$symbol, , drop = FALSE]) %*% s$coef)
mk <- function(time, event) { i <- if (shuffle) sample(length(time)) else seq_along(time); data.frame(time = time[i], event = event[i]) }

## Cohorts
sc <- readRDS("data/scanb_clin.rds"); sx <- readRDS("data/scanb_expr.rds")
tc <- readRDS("data/clin_tumor_survival.rds"); tx <- readRDS("data/expr_tumor_vst.rds")
sym2id <- setNames(gi$gene_id, gi$gene_name)
scanb_genes <- unique(c(sigs$A$symbol, sigs$B$symbol, hubs$symbol))
zs <- zscore(sx[intersect(scanb_genes, rownames(sx)), ])
tg <- unique(c(sigs$A$symbol, sigs$B$symbol, hubs$symbol)); tid <- sym2id[tg]; keep <- !is.na(tid) & tid %in% rownames(tx)
zt <- zscore(tx[tid[keep], ]); rownames(zt) <- tg[keep]

cohorts <- list(
  `SCAN-B` = list(y = mk(sc$OS_months, sc$OS_event), z = zs,
                  lp = list(A = rowMeans(cv$lpA), B = rowMeans(cv$lpB)), oof = TRUE),
  `TCGA-BRCA` = list(y = mk(tc$OS_months, tc$OS_event), z = zt,
                     lp = list(A = score(zt, sigs$A), B = score(zt, sigs$B)), oof = FALSE))
if (!shuffle && file.exists("data/metabric_prepared.rds")) {
  mb <- readRDS("data/metabric_prepared.rds")
  use <- function(s) { k <- s$symbol %in% rownames(mb$z); score(mb$z, s[k, ]) }
  cohorts$METABRIC <- list(y = mk(mb$clin$time, mb$clin$os), z = mb$z, lp = list(A = use(sigs$A), B = use(sigs$B)), oof = FALSE)
}

## Frozen development quantities from the final SCAN-B models: median cut-off and baseline hazard
final_lp <- lapply(sigs, function(s) score(zs, s))
cuts <- lapply(final_lp, median)
base <- lapply(final_lp, function(lp) {
  d <- cohorts$`SCAN-B`$y; d$lp <- lp
  b <- basehaz(coxph(Surv(time, event) ~ offset(lp), data = d), centered = FALSE)
  max_h <- function(t) { h <- b$hazard[b$time <= t]; if (length(h)) max(h) else 0 }
  max_h(CAL_T)
})

km_out <- roc_out <- cal_out <- list()
for (co in names(cohorts)) for (m in c("A", "B")) {
  ch <- cohorts[[co]]; d <- ch$y; d$lp <- ch$lp[[m]]
  cut_val <- if (ch$oof) median(d$lp) else cuts[[m]]
  d$group <- factor(ifelse(d$lp > cut_val, "High risk", "Low risk"), levels = c("Low risk", "High risk"))

  ## Kaplan-Meier and hazard ratio
  fit <- surv_fit(Surv(time, event) ~ group, data = d)
  cx  <- summary(coxph(Surv(time, event) ~ group, data = d))
  lr  <- survdiff(Surv(time, event) ~ group, data = d)
  km_out[[paste(co, m)]] <- data.frame(cohort = co, model = m, n_low = sum(d$group == "Low risk"), n_high = sum(d$group == "High risk"),
    events = sum(d$event), HR_high_vs_low = cx$conf.int[1, 1], HR_low = cx$conf.int[1, 3], HR_high = cx$conf.int[1, 4],
    logrank_p = 1 - pchisq(lr$chisq, 1))
  p <- ggsurvplot(fit, data = d, risk.table = TRUE, pval = TRUE, conf.int = FALSE, palette = c("#2166ac", "#b2182b"),
                  xlab = "Months", ylab = "Overall survival", title = sprintf("%s, Model %s", co, m),
                  legend.labs = c("Low risk", "High risk"), legend.title = "", xlim = c(0, 180), break.time.by = 30,
                  risk.table.height = 0.25, ggtheme = theme_bw(base_size = 11))
  png(file.path(figdir, sprintf("km_%s_model%s.png", gsub("[^A-Za-z]", "", co), m)), width = 1700 / 220, height = 1500 / 220, units = "in", res = 600)
  print(p); dev.off()

  ## Time-dependent ROC
  ok <- TIMES[TIMES < max(d$time[d$event == 1])]
  tr <- timeROC(T = d$time, delta = d$event, marker = d$lp, cause = 1, times = ok, iid = TRUE)
  ci <- confint(tr, level = 0.95)$CI_AUC / 100
  roc_out[[paste(co, m)]] <- data.frame(cohort = co, model = m, months = ok, at_risk = sapply(ok, function(t) sum(d$time >= t)),
    events_by_time = sapply(ok, function(t) sum(d$event == 1 & d$time <= t)),
    AUC = unname(tr$AUC[paste0("t=", ok)]), AUC_low = ci[, 1], AUC_high = ci[, 2], row.names = NULL)

  ## Calibration at the fixed horizon. SCAN-B: out-of-fold scores with their own baseline.
  ## External cohorts: frozen score with the SCAN-B baseline hazard.
  H0 <- if (ch$oof) {
    b <- basehaz(coxph(Surv(time, event) ~ offset(lp), data = d), centered = FALSE)
    h <- b$hazard[b$time <= CAL_T]; if (length(h)) max(h) else 0
  } else base[[m]]
  d$pred <- 1 - exp(-H0 * exp(d$lp))
  q <- unique(quantile(d$pred, probs = seq(0, 1, 0.2)))
  d$bin <- cut(d$pred, q, include.lowest = TRUE, labels = FALSE)
  obs <- do.call(rbind, lapply(sort(unique(d$bin)), function(b) {
    s <- summary(survfit(Surv(time, event) ~ 1, data = d[d$bin == b, ]), times = CAL_T, extend = TRUE)
    data.frame(bin = b, n = sum(d$bin == b), predicted = mean(d$pred[d$bin == b]), observed = 1 - s$surv, obs_low = 1 - s$upper, obs_high = 1 - s$lower)
  }))
  slope <- summary(coxph(Surv(time, event) ~ lp, data = d))
  cal_out[[paste(co, m)]] <- cbind(cohort = co, model = m, horizon_months = CAL_T, obs,
                                   calibration_slope = slope$coefficients[1, "coef"], slope_se = slope$coefficients[1, "se(coef)"])
  png(file.path(figdir, sprintf("calibration_%s_model%s.png", gsub("[^A-Za-z]", "", co), m)), width = 1300 / 220, height = 1300 / 220, units = "in", res = 600)
  lim <- c(0, max(obs$obs_high, obs$predicted, na.rm = TRUE) * 1.05)
  plot(obs$predicted, obs$observed, xlim = lim, ylim = lim, pch = 19, col = "#b2182b", xlab = sprintf("Predicted risk of death by %d months", CAL_T),
       ylab = "Observed risk (Kaplan-Meier)", main = sprintf("%s, Model %s", co, m))
  arrows(obs$predicted, obs$obs_low, obs$predicted, obs$obs_high, length = 0.03, angle = 90, code = 3, col = "#b2182b"); abline(0, 1, lty = 2, col = "grey40")
  dev.off()
}
km <- do.call(rbind, km_out); roc <- do.call(rbind, roc_out); cal <- do.call(rbind, cal_out)
rownames(km) <- rownames(roc) <- rownames(cal) <- NULL
suffix <- if (shuffle) "_SHUFFLED_TEST" else ""
write.csv(km, file.path(resdir, paste0("survival_km_summary", suffix, ".csv")), row.names = FALSE)
write.csv(roc, file.path(resdir, paste0("survival_time_roc", suffix, ".csv")), row.names = FALSE)
write.csv(cal, file.path(resdir, paste0("survival_calibration", suffix, ".csv")), row.names = FALSE)

## Hub genes: univariate Cox per SD of expression, within each cohort
hub_out <- do.call(rbind, lapply(names(cohorts), function(co) {
  ch <- cohorts[[co]]; present <- intersect(hubs$symbol, rownames(ch$z))
  out <- do.call(rbind, lapply(present, function(g) {
    s <- summary(coxph(Surv(ch$y$time, ch$y$event) ~ ch$z[g, ]))
    data.frame(cohort = co, symbol = g, module = hubs$module[match(g, hubs$symbol)], HR_per_SD = s$conf.int[1, 1],
               HR_low = s$conf.int[1, 3], HR_high = s$conf.int[1, 4], p = s$coefficients[1, "Pr(>|z|)"])
  }))
  out$FDR <- p.adjust(out$p, "BH"); out
}))
write.csv(hub_out, file.path(resdir, paste0("hub_gene_survival", suffix, ".csv")), row.names = FALSE)

cat("\nKaplan-Meier (high vs low risk, frozen SCAN-B cut-off for external cohorts):\n"); print(format(km, digits = 3), row.names = FALSE)
cat("\nTime-dependent AUC:\n"); print(format(roc, digits = 3), row.names = FALSE)
cat("\nCalibration slope at", CAL_T, "months:\n"); print(format(unique(cal[, c("cohort", "model", "calibration_slope", "slope_se")]), digits = 3), row.names = FALSE)
cat("\nHub genes with FDR < 0.05:\n"); print(format(hub_out[hub_out$FDR < 0.05, ], digits = 3), row.names = FALSE)
cat("\nDone.\n")
