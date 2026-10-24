#!/usr/bin/env Rscript
# =============================================================================
#  smoke_test.R -- dependency-free validation of the ChAMP2 core algorithms
#
#  Parses and sources every file in R/ and exercises the bundled algorithms
#  (KNN imputation, beta mixture EM, exact KPMT, kernel sums) without needing
#  ChAMPdata, limma or any other package.  This is the quickest way to check
#  that a source tree is compatible with an old or new R release:
#
#      Rscript tools/smoke_test.R
#
#  Exit status is 0 when every check passes.
# =============================================================================

root <- tryCatch({
    ofile <- sys.frame(1)$ofile
    if (is.null(ofile)) "." else normalizePath(file.path(dirname(ofile), ".."), mustWork = FALSE)
}, error = function(e) ".")
if (!dir.exists(file.path(root, "R"))) root <- normalizePath(".")
files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
if (!length(files)) stop("No R sources found; run this script from the package root.")

env <- new.env(parent = globalenv())
for (f in files) sys.source(f, envir = env)

checks <- 0L
failed <- 0L
check <- function(label, expr) {
    checks <<- checks + 1L
    ok <- tryCatch(isTRUE(expr), error = function(e) {
        cat("  FAIL ", label, ": ", conditionMessage(e), "\n", sep = "")
        FALSE
    })
    if (!ok) failed <<- failed + 1L else cat("  ok   ", label, "\n", sep = "")
    invisible(ok)
}

cat("ChAMP2 smoke test\n")
cat("  R          : ", R.version.string, "\n", sep = "")
cat("  sources    : ", length(files), " files\n", sep = "")
cat("  algorithms :\n")

## --- KNN imputation -------------------------------------------------------
set.seed(1)
x <- matrix(stats::rnorm(200), 20, 10,
            dimnames = list(paste0("p", 1:20), paste0("s", 1:10)))
x[sample(length(x), 20)] <- NA
y <- env$champ_knn_impute(x, k = 3, use_impute = FALSE, strict = FALSE)
check("KNN imputation fills every NA", !anyNA(y))
check("KNN imputation preserves dimnames", identical(dimnames(y), dimnames(x)))

## --- beta mixture ---------------------------------------------------------
set.seed(7)
b <- c(stats::rbeta(300, 2, 20), stats::rbeta(300, 20, 2), stats::rbeta(300, 5, 5))
w <- matrix(0, length(b), 3)
w[b < 0.2, 1] <- 1
w[b > 0.8, 3] <- 1
w[w[, 1] == 0 & w[, 3] == 0, 2] <- 1
fit <- env$champ_blc(matrix(b, ncol = 1), w, maxiter = 25, tol = 1e-6)
check("beta mixture separates low and high components", fit$mu[1, 1] < fit$mu[3, 1])
check("posterior weights sum to one", all(abs(rowSums(fit$w) - 1) < 1e-8))

## --- KPMT -----------------------------------------------------------------
pop <- stats::setNames(seq_len(500), paste0("g", seq_len(500)))
p_low <- env$champ_kpmt(pop, paste0("g", 1:20), tail = "lower")$p.value
p_high <- env$champ_kpmt(pop, paste0("g", 481:500), tail = "lower")$p.value
check("KPMT gives a tiny p-value for low ranks", p_low < 0.01)
check("KPMT gives a large p-value for high ranks", p_high > 0.9)

## --- kernel sums ----------------------------------------------------------
pos <- sort(sample.int(100000, 300))
k <- env$.champ_kernel_sums(pos, stats::runif(300), scale = 500)
check("kernel sums are positive", all(k[, "sky"] >= 0) && all(k[, "sk"] > 0))
check("kernel weights satisfy the Cauchy-Schwarz bound", all(k[, "skk"] <= k[, "sk"]^2 + 1e-8))

## --- annotation helpers ---------------------------------------------------
check("arraytype normalisation", identical(env$.champ_arraytype("450k"), "450K"))
check("genome mapping", identical(env$.champ_genome("EPICv2"), "hg38"))

cat("\n", if (failed == 0L) "ALL CHECKS PASSED" else paste0(failed, " CHECK(S) FAILED"),
    " (", checks, " run)\n", sep = "")
quit(status = if (failed == 0L) 0L else 1L)