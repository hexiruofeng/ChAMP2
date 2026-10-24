# =============================================================================
#  ChAMP 2.32 -- self-contained implementations of algorithms that used to be
#  imported from packages which are no longer reliably installable.
#
#  R 4.2.3 -- 4.6.2 spans Bioconductor 3.16 -- 3.24.  During that window the
#  following helper packages either disappeared, were archived or became hard
#  to install on current toolchains:
#
#    impute    (Bioconductor)  -> champ_knn_impute()   KNN imputation
#    RPMM      (CRAN archive)  -> champ_blc()          beta mixture (BMIQ)
#    kpmt      (CRAN archive)  -> champ_kpmt()         known-population median test
#    DMRcate   (still on Bioc) -> champ_dmrcate()      DMRcate compatible caller
#
#  Each function first uses the original package when it is installed (so
#  results are bit-identical to the legacy behaviour) and falls back to the
#  bundled implementation otherwise.  The fallbacks are written from the
#  published algorithms and are documented in ?champ_compat.
# =============================================================================

# ---------------------------------------------------------------------------
# numeric helpers
# ---------------------------------------------------------------------------

.champ_logsumexp <- function(x) {
    if (!length(x)) return(-Inf)
    m <- max(x)
    if (!is.finite(m)) return(m)
    m + log(sum(exp(x - m)))
}

#' @keywords internal
.champ_row_logsumexp <- function(m) {
    m <- as.matrix(m)
    mx <- apply(m, 1L, max)
    mx[!is.finite(mx)] <- 0
    mx + log(rowSums(exp(m - mx)))
}

# ===========================================================================
#  1. KNN imputation (replacement for impute::impute.knn)
# ===========================================================================

#' K-Nearest-Neighbour imputation
#'
#' Drop-in replacement for \code{impute::impute.knn}.  Rows are features
#' (probes) and columns are samples; missing entries of a feature are imputed
#' from the feature's values in the \code{k} most similar samples.
#'
#' The original \code{impute} package is used when available, so that existing
#' pipelines keep bit-identical results.  Otherwise a vectorised pure-R
#' implementation with the same arguments is used.
#'
#' @param data numeric matrix, rows = features, columns = samples.
#' @param k number of neighbours.
#' @param rowmax features with more than this fraction of missing values are
#'   filled with the column mean (as in \code{impute.knn}).
#' @param colmax samples with more than this fraction of missing values are
#'   rejected; for compatibility with the historical behaviour this is an
#'   error, unless \code{strict = FALSE}.
#' @param maxp maximum number of rows used to estimate sample distances
#'   (\code{impute.knn} splits the matrix instead; results are equivalent for
#'   the purpose of ChAMP).
#' @param use_impute use the \code{impute} package when installed (default).
#' @param strict error out on samples exceeding \code{colmax} (default TRUE,
#'   matching \code{impute.knn}); when FALSE such samples are mean-filled.
#' @param verbose print progress messages.
#'
#' @return the imputed matrix, with dimnames preserved.
#' @export
champ_knn_impute <- function(data, k = 5L, rowmax = 0.5, colmax = 0.8,
                             maxp = 1500L, use_impute = TRUE, strict = TRUE,
                             verbose = FALSE) {
    if (is.data.frame(data)) data <- as.matrix(data)
    if (!is.matrix(data) || !is.numeric(data)) {
        .champ_stop("'data' must be a numeric matrix.")
    }
    if (!anyNA(data)) return(data)
    dn <- dimnames(data)
    x <- matrix(as.numeric(data), nrow = nrow(data), ncol = ncol(data), dimnames = dn)

    if (isTRUE(use_impute) && .champ_has("impute")) {
        out <- .champ_try(impute::impute.knn(x, k = as.integer(k), rowmax = rowmax,
                                             colmax = colmax, maxp = as.integer(maxp))$data,
                          default = NULL)
        if (!is.null(out) && !anyNA(out)) return(out)
        if (is.null(out)) {
            .champ_msg("impute::impute.knn failed; using the ChAMP implementation instead.")
        }
    }

    n <- nrow(x)
    p <- ncol(x)
    if (p < 2L) {
        ## Single sample: use the feature mean.
        v <- x[, 1L]
        v[is.na(v)] <- mean(v, na.rm = TRUE)
        x[, 1L] <- v
        return(x)
    }
    k <- max(1L, min(as.integer(k), p - 1L))

    imiss <- is.na(x)
    col_nas <- colSums(imiss)
    bad_cols <- which(col_nas > colmax * n)
    if (length(bad_cols)) {
        msg <- paste0(sum(col_nas > colmax * n), " sample(s) have more than ",
                      round(colmax * 100), "% missing values (",
                      paste(colnames(x)[bad_cols], collapse = ", "), ").")
        if (isTRUE(strict)) {
            .champ_stop(msg, "\n  Use champ.impute(method = 'Delete') or lower ProbeCutoff first.")
        }
        .champ_warn(msg, " They will be filled with row means.")
    }

    ## rows with too many missing values: column-mean imputation
    row_na <- rowSums(imiss)
    too_many <- row_na > trunc(rowmax * p)
    if (any(too_many)) {
        .champ_msg("  ", sum(too_many), " feature(s) have more than ",
                   format(round(rowmax * 100, 1)), "% missing entries; column means are used.")
        cm <- colMeans(x[!too_many, , drop = FALSE], na.rm = TRUE)
        cm[is.na(cm)] <- 0.5
        for (j in seq_len(p)) x[too_many & imiss[, j], j] <- cm[j]
    }

    ## distance between samples, estimated from at most maxp features
    cand_rows <- which(row_na == 0L)
    if (!length(cand_rows)) cand_rows <- which(row_na <= trunc(rowmax * p))
    if (length(cand_rows) > maxp) {
        set.seed(1234567)
        cand_rows <- sort(sample(cand_rows, maxp))
    }
    sub <- x[cand_rows, , drop = FALSE]
    ## remaining NAs in the distance block are filled with the column mean of the block
    if (anyNA(sub)) {
        cm2 <- colMeans(sub, na.rm = TRUE)
        cm2[is.na(cm2)] <- 0.5
        sub[is.na(sub)] <- rep(cm2, each = nrow(sub))[is.na(sub)]
    }
    d <- as.matrix(stats::dist(t(sub)))
    d[!is.finite(d)] <- max(d[is.finite(d)], na.rm = TRUE)

    rows_na <- which(row_na > 0L & !too_many)
    if (length(rows_na)) {
        ## group rows by missingness pattern so the neighbour search is shared
        patt <- apply(imiss[rows_na, , drop = FALSE], 1L,
                      function(z) paste(which(z), collapse = ","))
        groups <- split(rows_na, patt)
        if (isTRUE(verbose)) .champ_msg("  ", length(groups), " missingness pattern(s).")
        for (g in groups) {
            miss <- which(imiss[g[1L], ])
            obs  <- setdiff(seq_len(p), miss)
            if (!length(obs)) {
                ## nothing observed for this pattern: column mean
                cm3 <- colMeans(x[-g, , drop = FALSE], na.rm = TRUE)
                cm3[is.na(cm3)] <- 0.5
                for (j in miss) x[g, j] <- cm3[j]
                next
            }
            donors <- obs[col_nas[obs] <= colmax * n]
            if (!length(donors)) donors <- obs
            for (j in miss) {
                dd <- d[j, donors]
                nn <- donors[order(dd)[seq_len(min(k, length(donors)))]]
                x[g, j] <- rowMeans(x[g, nn, drop = FALSE], na.rm = TRUE)
            }
        }
    }
    attr(x, "imputation") <- "champ_knn_impute"
    x
}

# ===========================================================================
#  2. Beta mixture EM (replacement for RPMM::blc, used by BMIQ)
# ===========================================================================

#' Maximum likelihood beta parameters for one component
#' @keywords internal
.champ_beta_est <- function(y, w, weights) {
    yobs <- !is.na(y)
    if (sum(yobs) <= 1L) return(c(1, 1))
    y <- y[yobs]; w <- w[yobs]; weights <- weights[yobs]
    N <- sum(weights * w)
    if (!is.finite(N) || N <= 0) return(c(1, 1))
    m <- sum(weights * w * y) / N
    v <- sum(weights * w * y * y) / N - m * m
    if (!is.finite(v) || v <= 1e-10) v <- 1e-10
    m <- min(max(m, 1e-6), 1 - 1e-6)
    conc <- max(1e-6, m * (1 - m) / v - 1)
    logab <- log(c(m, 1 - m)) + log(conc)
    if (sum(yobs) == 2L) return(exp(logab))
    nll <- function(la) {
        ab <- exp(la)
        -sum(w * weights * stats::dbeta(y, ab[1L], ab[2L], log = TRUE), na.rm = TRUE)
    }
    opt <- try(stats::optim(logab, nll, method = "BFGS"), silent = TRUE)
    if (inherits(opt, "try-error") || !is.finite(opt$value)) return(exp(logab))
    par <- exp(opt$par)
    if (any(!is.finite(par)) || any(par <= 0)) return(exp(logab))
    par
}

#' Beta latent class model (EM)
#'
#' Re-implementation of the beta mixture EM used by BMIQ.  Returns the same
#' components as \code{RPMM::blc}: \code{a}, \code{b}, \code{eta}, \code{mu},
#' \code{w} and the log-likelihood.
#'
#' @param Y numeric matrix (n x J), typically a single column of beta values.
#' @param w numeric matrix (n x K) of initial posterior class memberships.
#' @param maxiter maximum EM iterations.
#' @param tol convergence tolerance on mu.
#' @param weights case weights.
#' @param verbose print the convergence criterion at each iteration.
#'
#' @return list with a, b (K x J), eta (K), mu (K x J), w (n x K), llike.
#' @export
champ_blc <- function(Y, w, maxiter = 25L, tol = 1e-6, weights = NULL,
                      verbose = FALSE) {
    Y <- as.matrix(Y)
    w <- as.matrix(w)
    if (nrow(Y) != nrow(w)) .champ_stop("Dimensions of 'w' and 'Y' do not agree.")
    n <- nrow(Y); J <- ncol(Y); K <- ncol(w)
    if (is.null(weights)) weights <- rep(1, n)

    Ymn <- min(Y[Y > 0], na.rm = TRUE)
    Ymx <- max(Y[Y < 1], na.rm = TRUE)
    if (!is.finite(Ymn)) Ymn <- 1e-6
    if (!is.finite(Ymx)) Ymx <- 1 - 1e-6
    Y <- pmax(Y, Ymn / 2)
    Y <- pmin(Y, 1 - (1 - Ymx) / 2)
    Yobs <- !is.na(Y)

    mu <- a <- b <- matrix(Inf, K, J)
    eta <- rep(1 / K, K)
    crit <- Inf
    for (iter in seq_len(maxiter)) {
        eta <- colSums(weights * w) / sum(weights)
        mu0 <- mu
        for (k in seq_len(K)) {
            for (j in seq_len(J)) {
                ab <- .champ_beta_est(Y[, j], w[, k], weights)
                a[k, j] <- ab[1L]
                b[k, j] <- ab[2L]
                mu[k, j] <- ab[1L] / (ab[1L] + ab[2L])
            }
        }
        ww <- array(0, dim = c(n, J, K))
        for (k in seq_len(K)) {
            for (j in seq_len(J)) {
                ww[Yobs[, j], j, k] <- stats::dbeta(Y[Yobs[, j], j], a[k, j], b[k, j],
                                                     log = TRUE)
            }
        }
        ll <- apply(ww, c(1L, 3L), sum, na.rm = TRUE)
        mx <- apply(ll, 1L, max)
        llc <- ll - mx
        num <- t(eta * t(exp(llc)))
        like <- rowSums(num)
        w <- num / like
        llike <- weights * (log(like) + mx)
        crit <- max(abs(mu - mu0))
        if (isTRUE(verbose)) .champ_msg("    blc iter ", iter, " crit ", signif(crit, 3))
        if (is.finite(crit) && crit < tol) break
    }
    list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike, na.rm = TRUE))
}

# ===========================================================================
#  3. Known population median test (replacement for kpmt::kpmt)
# ===========================================================================

#' Known Population Median Test
#'
#' Exact finite-population test for whether the minimal median of an observed
#' subset is shifted relative to the whole population.  This is a
#' re-implementation of \code{kpmt::kpmt} (CRAN, archived 2025) and returns a
#' data frame with the same columns, so existing code that indexes the fourth
#' column (\code{p.value}) keeps working.
#'
#' The p-value is the exact combinatorial probability
#' \deqn{P(\tilde{X} \le w) = \frac{s!(N-s)!}{N!(m-1)!(s-m)!}
#'       \sum_{x \in U} \frac{(x-1)!(N-x)!}{(x-m)!(N-x-s+m)!}}
#' with \eqn{m = \lceil s/2 \rceil} and \eqn{U} the population ranks that are
#' below the observed minimal median and still leave room for the sample.
#'
#' @param pop named numeric vector (or one-column matrix) of population values.
#' @param obs character vector of population member names that were observed.
#' @param size sample size (only needed when \code{med} is supplied).
#' @param med pre-computed minimal median.
#' @param tail one of "lower", "upper", "two-sided".
#'
#' @return data.frame with columns feature, sample, logp, p.value,
#'   median.sample, median.all, median.other.
#' @export
champ_kpmt <- function(pop, obs = NULL, med = NULL, size = NULL,
                       tail = "lower") {
    if (is.matrix(pop) || is.data.frame(pop)) {
        df <- as.data.frame(pop)
        nms <- rownames(df)
        pop <- df[[1L]]
    } else {
        nms <- names(pop)
    }
    pop <- as.numeric(pop)
    if (is.null(nms) || length(nms) != length(pop) || anyNA(nms)) {
        nms <- as.character(seq_along(pop))
    }
    names(pop) <- nms
    N <- length(pop)
    if (is.null(obs) && is.null(med)) .champ_stop("Provide either 'obs' or 'med'.")

    if (!is.null(obs)) {
        obs <- as.character(obs)
        if (!all(obs %in% nms)) {
            .champ_stop("Observations not found in the population: ",
                        paste(utils::head(setdiff(obs, nms), 5L), collapse = ", "))
        }
        s <- length(obs)
        vals <- pop[match(obs, nms)]
        med_obs <- sort(vals)[ceiling(s / 2)]
    } else {
        med_obs <- as.numeric(med)[1L]
        if (is.null(size)) .champ_stop("'size' is required when 'med' is given.")
        s <- as.integer(size[1L])
    }
    if (s < 2L) .champ_stop("The sample must contain at least two members.")
    m <- ceiling(s / 2)

    p_from <- function(values, w, lower = TRUE) {
        v <- sort(values)
        if (!lower) { v <- -rev(v); w <- -w }
        x <- seq_len(N)
        logsum <- lgamma(x) + lgamma(N - x + 1) -
            lgamma(x - m + 1) - lgamma(N - x - s + m + 1)
        U <- which(v <= w & x >= m & (N - x) >= (s - m))
        if (!length(U)) return(-Inf)
        const <- lgamma(s + 1) + lgamma(N - s + 1) - lgamma(N + 1) -
            lgamma(m) - lgamma(s - m + 1)
        const + .champ_logsumexp(logsum[U])
    }

    tl <- tolower(tail[1L])
    logp_low <- p_from(pop, med_obs, lower = TRUE)
    logp_up  <- p_from(pop, med_obs, lower = FALSE)
    logp <- switch(tl,
        "lower"     = logp_low,
        "upper"     = logp_up,
        "two-sided" = .champ_logsumexp(c(logp_low, logp_up)),
        .champ_stop("'tail' must be one of lower, upper, two-sided."))
    logp <- min(logp, 0)

    data.frame(
        feature       = "feature1",
        sample        = "sample1",
        logp          = logp,
        p.value       = exp(logp),
        median.sample = med_obs,
        median.all    = sort(pop)[ceiling(N / 2)],
        median.other  = NA_real_,
        stringsAsFactors = FALSE)
}

# ===========================================================================
#  4. DMRcate-compatible differential methylation regions
# ===========================================================================

#' Weighted sum of values per bin, without tabulate(weights = )
#'
#' tabulate() lost its weights argument in R 4.6, so the binning is done with
#' rowsum(), which is equally fast and available on every supported release.
#' @keywords internal
.champ_weighted_sum <- function(bins, weights, n) {
    out <- numeric(n)
    if (!length(bins)) return(out)
    s <- rowsum(as.numeric(weights), bins, reorder = TRUE)
    idx <- as.integer(rownames(s))
    ok <- !is.na(idx) & idx >= 1L & idx <= n
    out[idx[ok]] <- s[ok, 1L]
    out
}

#' Kernel sums used by the DMRcate test
#' @keywords internal
.champ_kernel_sums <- function(pos, y, scale, support = 5) {
    n <- length(pos)
    if (n < 2L) return(cbind(sky = y, sk = rep(1, n), skk = rep(1, n)))
    ord <- order(pos)
    pos_o <- pos[ord]; y_o <- y[ord]
    maxdiff <- support * scale
    j <- findInterval(pos_o + maxdiff, pos_o)
    d <- j - seq_len(n)
    keep <- d > 0L
    sky <- y_o; sk <- rep(1, n); skk <- rep(1, n)
    if (any(keep)) {
        ii <- rep(seq_len(n)[keep], d[keep])
        jj <- ii + sequence(d[keep])
        dx <- pos_o[jj] - pos_o[ii]
        k <- exp(-(dx / scale)^2 / 2)
        sky <- sky + .champ_weighted_sum(ii, k * y_o[jj], n) +
            .champ_weighted_sum(jj, k * y_o[ii], n)
        sk <- sk + .champ_weighted_sum(ii, k, n) + .champ_weighted_sum(jj, k, n)
        skk <- skk + .champ_weighted_sum(ii, k^2, n) +
            .champ_weighted_sum(jj, k^2, n)
    }
    out <- cbind(sky = sky, sk = sk, skk = skk)
    out[order(ord), , drop = FALSE]
}

#' Satterthwaite chi-square kernel test
#' @keywords internal
.champ_kernel_test <- function(pos, X2, lambda, df = 1) {
    j <- .champ_kernel_sums(pos, X2, scale = lambda)
    sk <- j[, "sk"]; skk <- j[, "skk"]; sky <- j[, "sky"]
    a <- df * skk / sk
    b <- sk * sk / skk
    stats::pchisq(sky / a, df = b, lower.tail = FALSE)
}

#' DMRcate-compatible DMR calling
#'
#' A self-contained implementation of the DMRcate algorithm
#' (Peters et al., Epigenetics & Chromatin 2015): probes are annotated with
#' moderated t-statistics, a Gaussian kernel with bandwidth
#' \code{lambda / C} is used to compute local chi-square statistics, and
#' regions of consecutive significant probes are summarised with Stouffer,
#' harmonic-mean and Fisher combination of p-values.
#'
#' When the \code{DMRcate} package is installed, \code{champ.DMR()} uses it by
#' default; this function is the fallback, and produces a data frame with the
#' same columns as \code{as.data.frame(DMRcate::extractRanges(...))}.
#'
#' @param beta beta value matrix (probes x samples).
#' @param pheno phenotype vector with exactly two groups (or a numeric
#'   covariate).
#' @param arraytype one of "450K", "EPICv1", "EPICv2".
#' @param fdr probe-level FDR used to define significant probes.
#' @param min.cpgs minimum number of probes per DMR.
#' @param lambda kernel bandwidth in base pairs.
#' @param C scaling of the bandwidth (\code{sigma = lambda / C}).
#' @param pcutoff "fdr" or a numeric probe-level p-value cutoff.
#' @param betacutoff optional absolute mean beta difference filter.
#' @param verbose print progress.
#'
#' @return data.frame with columns seqnames, start, end, width, strand,
#'   no.cpgs, min_smoothed_fdr, Stouffer, HMFDR, Fisher, maxdiff, meandiff,
#'   overlapping.genes.
#' @export
champ_dmrcate <- function(beta, pheno, arraytype = "450K", fdr = 0.05,
                          min.cpgs = 7L, lambda = 1000, C = 2,
                          pcutoff = "fdr", betacutoff = NULL, verbose = TRUE) {
    .champ_require("limma", "champ_dmrcate()")
    beta <- .champ_as_beta(beta)
    at <- .champ_arraytype(arraytype)
    pheno <- .champ_pheno(pheno, ncol(beta))
    if (length(unique(pheno)) != 2L) {
        .champ_stop("champ_dmrcate() needs exactly two groups in 'pheno'.")
    }

    pf <- .champ_probe_features(at)
    keep <- rownames(beta) %in% rownames(pf)
    if (!any(keep)) .champ_stop("None of the probes in 'beta' are present in the annotation.")
    beta <- beta[keep, , drop = FALSE]
    anno <- pf[rownames(beta), , drop = FALSE]
    ok <- !is.na(anno$MAPINFO) & !is.na(anno$CHR) & anno$CHR %in% c(as.character(1:26), "X", "Y")
    beta <- beta[ok, , drop = FALSE]
    anno <- anno[ok, , drop = FALSE]
    if (nrow(beta) < min.cpgs) .champ_stop("Too few annotated probes to call DMRs.")

    M <- .champ_logit2(beta)
    design <- stats::model.matrix(~ factor(pheno))
    coef <- ncol(design)
    fit <- limma::eBayes(limma::lmFit(M, design))
    tt <- limma::topTable(fit, coef = coef, number = nrow(M), sort.by = "none")
    bfit <- limma::eBayes(limma::lmFit(beta, design))
    btt <- limma::topTable(bfit, coef = coef, number = nrow(beta), sort.by = "none")
    m <- match(rownames(tt), rownames(btt))

    object <- data.frame(
        ID = rownames(tt),
        weights = abs(tt$t),
        CHR = anno$CHR,
        pos = as.numeric(anno$MAPINFO),
        rawpval = tt$P.Value,
        diff = btt$logFC[m],
        indfdr = tt$adj.P.Val,
        is.sig = tt$adj.P.Val < fdr,
        stringsAsFactors = FALSE)
    object <- object[order(object$CHR, object$pos), , drop = FALSE]
    tab <- table(object$CHR)
    object <- object[!object$CHR %in% names(tab)[tab == 1L], , drop = FALSE]
    if (!nrow(object)) .champ_stop("No DMR could be called (not enough probes per chromosome).")

    ## --- per chromosome kernel test -------------------------------------
    sigma <- lambda / C
    object$raw <- NA_real_
    for (chr in unique(object$CHR)) {
        idx <- which(object$CHR == chr)
        if (length(idx) < 2L) next
        if (isTRUE(verbose)) .champ_msg("    kernel test: chr", chr, " (", length(idx), " probes)")
        object$raw[idx] <- .champ_kernel_test(object$pos[idx], object$weights[idx]^2,
                                              lambda = sigma, df = 1)
    }
    object <- object[!is.na(object$raw), , drop = FALSE]
    object$fdr <- stats::p.adjust(object$raw, method = "BH")

    if (identical(pcutoff, "fdr")) {
        nsig <- sum(object$is.sig)
        if (nsig == 0L) {
            .champ_stop("No probe passed the FDR threshold in cpg.annotate(); ",
                        "try increasing 'fdr' or set 'pcutoff' explicitly.")
        }
        pcutoff <- sort(object$fdr)[nsig]
    }
    object$sig <- object$fdr <= pcutoff
    if (sum(object$sig) < 2L) .champ_stop("Fewer than two significant probes; no DMR can be formed.")

    ## --- segmentation ---------------------------------------------------
    sig_idx <- which(object$sig)
    pos_k <- object$pos[sig_idx]
    chr_k <- object$CHR[sig_idx]
    jump <- (chr_k[-1L] != chr_k[-length(chr_k)]) | (diff(pos_k) > lambda)
    breaks <- c(1L, which(jump) + 1L, length(sig_idx) + 1L)
    regions <- lapply(seq_len(length(breaks) - 1L), function(i) {
        sig_idx[seq(breaks[i], breaks[i + 1L] - 1L)]
    })
    if (!length(regions)) .champ_stop("No DMR could be formed.")

    fn_stouffer <- function(x) stats::pnorm(sum(stats::qnorm(x)) / sqrt(length(x)))
    fn_hm <- function(x) 1 / mean(1 / x)
    fn_fisher <- function(x) stats::pchisq(-2 * sum(log(x)), df = length(x) * 2, lower.tail = FALSE)
    fn_max <- function(x) x[which.max(abs(x))]

    res <- lapply(regions, function(idx) {
        data.frame(
            seqnames = object$CHR[idx[1L]],
            start = min(object$pos[idx]),
            end = max(object$pos[idx]),
            no.cpgs = length(idx),
            min_smoothed_fdr = min(object$indfdr[idx]),
            Stouffer = fn_stouffer(object$rawpval[idx]),
            HMFDR = fn_hm(object$rawpval[idx]),
            Fisher = fn_fisher(object$rawpval[idx]),
            maxdiff = fn_max(object$diff[idx]),
            meandiff = mean(object$diff[idx]),
            stringsAsFactors = FALSE)
    })
    res <- do.call(rbind, res)
    res$Stouffer <- stats::p.adjust(res$Stouffer, method = "BH")
    res$Fisher <- stats::p.adjust(res$Fisher, method = "BH")
    res$HMFDR <- stats::p.adjust(res$HMFDR, method = "BH")
    res <- res[res$no.cpgs >= min.cpgs, , drop = FALSE]
    if (!is.null(betacutoff)) res <- res[abs(res$meandiff) > betacutoff, , drop = FALSE]
    if (!nrow(res)) .champ_stop("No DMR passed the min.cpgs/betacutoff filters.")
    res <- res[order(res$min_smoothed_fdr, -res$no.cpgs), , drop = FALSE]

    res$width <- res$end - res$start + 1L
    res$strand <- "*"
    res$overlapping.genes <- .champ_genes_in_regions(res, pf, at)
    res <- res[, c("seqnames", "start", "end", "width", "strand", "no.cpgs",
                   "min_smoothed_fdr", "Stouffer", "HMFDR", "Fisher",
                   "maxdiff", "meandiff", "overlapping.genes")]
    rownames(res) <- paste0("DMR_", seq_len(nrow(res)))
    res
}

#' Genes overlapping DMR intervals, from the ChAMP annotation.
#' @keywords internal
.champ_genes_in_regions <- function(dmr, pf, arraytype = "450K") {
    out <- rep(NA_character_, nrow(dmr))
    for (i in seq_len(nrow(dmr))) {
        sel <- pf$CHR == dmr$seqnames[i] & pf$MAPINFO >= dmr$start[i] &
            pf$MAPINFO <= dmr$end[i] & nzchar(pf$gene)
        if (!any(sel, na.rm = TRUE)) next
        g <- unlist(strsplit(pf$gene[which(sel)], ";", fixed = TRUE), use.names = FALSE)
        g <- unique(trimws(g))
        g <- g[nzchar(g)]
        if (length(g)) out[i] <- paste(g, collapse = ", ")
    }
    out
}