# =============================================================================
#  DoPBC -- Peak-Based Correction (Dedeurwaerder et al., 2013)
# =============================================================================

#' Peak-based correction of type II probe bias
#'
#' @param beta.m beta value matrix (probes x samples).
#' @param design.v probe design: 1 = type I, 2 = type II.
#' @param verbose print progress.
#'
#' @return corrected beta value matrix.
#' @export
DoPBC <- function(beta.m, design.v, verbose = TRUE) {
    beta.m <- as.matrix(beta.m)
    if (length(design.v) != nrow(beta.m)) {
        .champ_stop("design.v and beta.m have incompatible dimensions.")
    }
    eps <- 1e-6
    b <- pmin(pmax(beta.m, eps), 1 - eps)
    mval.m <- log2(b / (1 - b))
    type1.idx <- which(design.v == 1)
    type2.idx <- which(design.v == 2)
    mvalT.m <- mval.m
    for (s in seq_len(ncol(beta.m))) {
        neg.idx <- which(mval.m[, s] < 0)
        pos.idx <- which(mval.m[, s] > 0)
        neg1.idx <- intersect(neg.idx, type1.idx)
        neg2.idx <- intersect(neg.idx, type2.idx)
        pos1.idx <- intersect(pos.idx, type1.idx)
        pos2.idx <- intersect(pos.idx, type2.idx)

        peak <- function(idx) {
            if (length(idx) < 2L) return(NA_real_)
            d.o <- stats::density(mval.m[idx, s], kernel = "gaussian", bw = 0.5)
            abs(d.o$x[which.max(d.o$y)])
        }
        peakU2 <- peak(neg2.idx)
        peakM2 <- peak(pos2.idx)
        peakU1 <- peak(neg1.idx)
        peakM1 <- peak(pos1.idx)

        if (length(neg2.idx) && is.finite(peakU2) && peakU2 > 0 && is.finite(peakU1)) {
            mvalT.m[neg2.idx, s] <- (mval.m[neg2.idx, s] / peakU2) * peakU1
        }
        if (length(pos2.idx) && is.finite(peakM2) && peakM2 > 0 && is.finite(peakM1)) {
            mvalT.m[pos2.idx, s] <- (mval.m[pos2.idx, s] / peakM2) * peakM1
        }
        if (isTRUE(verbose)) .champ_msg("Done for sample ", s)
    }
    betaT.m <- 2^mvalT.m / (2^mvalT.m + 1)
    dimnames(betaT.m) <- dimnames(beta.m)
    betaT.m
}
