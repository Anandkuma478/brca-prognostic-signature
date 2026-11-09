suppressPackageStartupMessages({ library(survival); library(glmnet); library(doParallel) })
registerDoParallel(cores = 8)
set.seed(2028)
dir.create("results", showWarnings = FALSE)

sc <- readRDS("data/scanb_clin.rds"); sx <- readRDS("data/scanb_expr.rds")
cand <- readRDS("data/signature_scanb_model.rds")$candidates
sigB <- read.csv("results/signature_scanb_sparse30.csv")
Xs <- scale(t(sx[cand$symbol, ]))
ys <- Surv(sc$OS_months, sc$OS_event)
n <- nrow(Xs); B <- 100; MAXG <- 30

idx_list <- lapply(seq_len(B), function(b) {
  unlist(lapply(c(0, 1), function(e) { i <- which(sc$OS_event == e); sample(i, length(i), replace = TRUE) }))
})

fit_one <- function(idx) {
  fit <- suppressWarnings(glmnet(Xs[idx, ], ys[idx], family = "cox", alpha = 1, dfmax = MAXG + 10))
  k <- max(which(fit$df <= MAXG))
  b <- fit$beta[, k]
  b[b != 0]
}

t0 <- Sys.time()
first <- fit_one(idx_list[[1]])
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Bootstrap 1 took %.0f s and selected %d genes. Estimated total (8 cores): ~%.0f min\n",
            dt, length(first), dt * B / 8 / 60))

res <- foreach(b = 2:B, .packages = c("glmnet", "survival")) %dopar% fit_one(idx_list[[b]])
res <- c(list(first), res)
cat("All", B, "bootstraps done\n")

genes <- cand$symbol
sel  <- sapply(genes, function(g) sum(sapply(res, function(r) g %in% names(r))))
mcoef <- sapply(genes, function(g) { v <- unlist(lapply(res, function(r) if (g %in% names(r)) r[[g]])); if (length(v)) mean(v) else NA })
sgn   <- sapply(genes, function(g) { v <- unlist(lapply(res, function(r) if (g %in% names(r)) r[[g]])); if (length(v)) max(mean(v > 0), mean(v < 0)) else NA })
out <- data.frame(symbol = genes, selected = sel, freq = sel / B, mean_coef = mcoef, sign_consistency = sgn)
out <- out[order(-out$freq), ]
write.csv(out, "results/stability_selection.csv", row.names = FALSE)

cat("\nGenes selected at least once:", sum(out$freq > 0), "of", nrow(out), "\n")
cat("Selection frequency >= 0.5:", sum(out$freq >= 0.5), "| >= 0.8:", sum(out$freq >= 0.8), "\n")
cat("Median genes selected per bootstrap:", median(sapply(res, length)), "\n")
cat("\nTop 30 genes by selection frequency:\n"); print(head(out, 30), row.names = FALSE)
cat("\nThe", nrow(sigB), "genes of the model-B final signature, with their bootstrap frequency:\n")
print(out[match(sigB$symbol, out$symbol), ], row.names = FALSE)
cat("\nOf the", nrow(sigB), "final-signature genes, frequency >= 0.5:", sum(out$freq[match(sigB$symbol, out$symbol)] >= 0.5), "\n")
cat("\nSaved: results/stability_selection.csv\n")
