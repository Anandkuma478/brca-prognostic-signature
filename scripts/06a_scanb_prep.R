if (!requireNamespace("data.table", quietly = TRUE)) stop("Package missing: data.table")
library(data.table)
dir.create("data", showWarnings = FALSE)
gdir <- "data/geo"
expr_file <- file.path(gdir, "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")

read_series <- function(f) {
  L <- readLines(gzfile(f), warn = FALSE)
  L <- L[startsWith(L, "!Sample_")]
  sp <- strsplit(L, "\t", fixed = TRUE)
  clean <- function(x) gsub('^"|"$', "", x)
  title <- clean(sp[[which(startsWith(L, "!Sample_title"))[1]]][-1])
  gsm <- clean(sp[[which(startsWith(L, "!Sample_geo_accession"))[1]]][-1])
  out <- data.frame(gsm = gsm, title = title, stringsAsFactors = FALSE)
  for (i in which(startsWith(L, "!Sample_characteristics_ch1"))) {
    v <- clean(sp[[i]][-1])
    key <- sub(":.*$", "", v[which(nzchar(v))[1]])
    val <- sub("^[^:]*: ?", "", v)
    if (length(val) == nrow(out) && !(key %in% names(out))) out[[key]] <- val
  }
  out
}
a <- read_series(file.path(gdir, "GSE96058-GPL11154_series_matrix.txt.gz"))
b <- read_series(file.path(gdir, "GSE96058-GPL18573_series_matrix.txt.gz"))
cn <- intersect(names(a), names(b))
clin <- rbind(a[, cn], b[, cn])
cat("Samples in series matrices:", nrow(clin), "\n")
cat("Example titles:", head(clin$title, 5), "\n")

need <- c("overall survival days", "overall survival event", "pam50 subtype", "age at diagnosis",
          "tumor size", "lymph node status", "er status", "her2 status", "nhg",
          "endocrine treated", "chemo treated")
stopifnot(all(need %in% names(clin)))

hdr <- names(fread(cmd = paste("zcat", shQuote(expr_file)), nrows = 0))
cat("Expression columns:", length(hdr) - 1, " examples:", head(hdr[-1], 5), "\n")

is_repl <- grepl("repl", clin$title, ignore.case = TRUE)
cat("Replicate samples (dropped):", sum(is_repl), "\n")
clin <- clin[!is_repl, ]
clin$col <- ifelse(clin$title %in% hdr, clin$title, NA)
cat("Primary samples matched to expression columns:", sum(!is.na(clin$col)), "of", nrow(clin), "\n")
if (mean(!is.na(clin$col)) < 0.9) stop("Titles do not match expression columns: paste the 'Example titles' and 'Expression columns' lines")
clin <- clin[!is.na(clin$col), ]

num <- function(x) suppressWarnings(as.numeric(x))
clin$OS_days   <- num(clin[["overall survival days"]])
clin$OS_event  <- num(clin[["overall survival event"]])
clin$OS_months <- clin$OS_days / 30.44
clin$age       <- num(clin[["age at diagnosis"]])
clin$tumor_size_mm <- num(clin[["tumor size"]])
clin$node_pos  <- ifelse(clin[["lymph node status"]] == "NodePositive", 1,
                  ifelse(clin[["lymph node status"]] == "NodeNegative", 0, NA))
clin$er   <- num(clin[["er status"]])
clin$her2 <- num(clin[["her2 status"]])
clin$nhg  <- factor(ifelse(clin[["nhg"]] %in% c("G1", "G2", "G3"), clin[["nhg"]], NA))
clin$pam50 <- factor(ifelse(clin[["pam50 subtype"]] %in% c("Basal", "Her2", "LumA", "LumB", "Normal"),
                            clin[["pam50 subtype"]], NA))
clin$endocrine <- num(clin[["endocrine treated"]])
clin$chemo     <- num(clin[["chemo treated"]])

bad <- is.na(clin$OS_days) | is.na(clin$OS_event) | clin$OS_days <= 0
cat("Removed (missing or non-positive survival):", sum(bad), "\n")
clin <- clin[!bad, ]
rownames(clin) <- clin$col

cols <- c(1, match(clin$col, hdr))
dt <- fread(cmd = paste("zcat", shQuote(expr_file)), select = cols)
m <- as.matrix(dt[, -1]); rownames(m) <- dt[[1]]; colnames(m) <- hdr[cols[-1]]
if (anyDuplicated(rownames(m))) {
  m <- m[order(rowMeans(m), decreasing = TRUE), ]
  m <- m[!duplicated(rownames(m)), ]
}
m <- m[, rownames(clin)]
cat("\nExpression matrix (genes x samples):", dim(m), "\n")
cat("Value range:", round(range(m), 2), "\n")

cat("\nFinal SCAN-B cohort:", nrow(clin), "patients\n")
cat("Deaths:", sum(clin$OS_event), " Censored:", sum(clin$OS_event == 0), "\n")
cat("Median observed time (months):", round(median(clin$OS_months), 1), "\n")
cat("PAM50:\n"); print(table(clin$pam50, useNA = "ifany"))
cat("ER status (1=pos):\n"); print(table(clin$er, useNA = "ifany"))
cat("Node positive:\n"); print(table(clin$node_pos, useNA = "ifany"))
cat("Grade:\n"); print(table(clin$nhg, useNA = "ifany"))
cat("Chemo / endocrine treated:\n"); print(table(clin$chemo, clin$endocrine, useNA = "ifany"))
cat("Age summary:\n"); print(summary(clin$age))
cat("Genes in TCGA vs SCAN-B symbol overlap:\n")
gi <- tryCatch(as.data.frame(readRDS("data/gene_info.rds")), error = function(e) NULL)
if (!is.null(gi)) {
  nm <- intersect(c("gene_name", "external_gene_name", "symbol"), colnames(gi))[1]
  cat(" TCGA gene symbols found in SCAN-B:", sum(unique(gi[[nm]]) %in% rownames(m)), "of", length(unique(gi[[nm]])), "\n")
}
saveRDS(clin, "data/scanb_clin.rds"); saveRDS(m, "data/scanb_expr.rds")
cat("\nSaved: data/scanb_clin.rds, data/scanb_expr.rds\n")
