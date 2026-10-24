# =============================================================================
#  champ.refbase -- reference based cell proportion estimation
# =============================================================================

#' Estimate and correct blood cell proportions
#'
#' Uses a reference-based (constrained least squares) deconvolution with the
#' Reinius et al. (450K) or Accomando et al. (27K) reference profiles and
#' regresses out the estimated proportions.
#'
#' @param beta beta value matrix.
#' @param arraytype "450K", "EPICv1" (both use the 450K reference) or "27K".
#' @param verbose print progress.
#'
#' @return list(CorrectedBeta, CellFraction).
#' @export
champ.refbase <- function(beta = NULL, arraytype = "450K", verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<< ChAMP.REFBASE START >>>]")
    .champ_msg("-----------------------------")
    .champ_require("quadprog", "champ.refbase()")
    beta <- .champ_as_beta(beta)
    at <- .champ_arraytype(arraytype)
    if (at %in% c("EPICv1", "EPICv2")) {
        .champ_msg("Using the 450K reference profiles for ", at, ".")
    }

    projectWBC <- function(Y, coefWBC, nonnegative = TRUE, lessThanOne = FALSE) {
        Xmat <- as.matrix(coefWBC)
        nCol <- ncol(Xmat)
        nSubj <- ncol(Y)
        mixCoef <- matrix(0, nSubj, nCol, dimnames = list(colnames(Y), colnames(Xmat)))
        for (i in seq_len(nSubj)) {
            obs <- which(!is.na(Y[, i]))
            if (length(obs) < nCol) next
            Dmat <- crossprod(Xmat[obs, , drop = FALSE])
            dvec <- crossprod(Xmat[obs, , drop = FALSE], Y[obs, i])
            if (isTRUE(nonnegative)) {
                if (isTRUE(lessThanOne)) {
                    Amat <- cbind(rep(-1, nCol), diag(nCol))
                    b0vec <- c(-1, rep(0, nCol))
                } else {
                    Amat <- diag(nCol)
                    b0vec <- rep(0, nCol)
                }
                sol <- .champ_try(quadprog::solve.QP(Dmat, dvec, Amat, b0vec)$sol,
                                  default = NULL)
                if (is.null(sol)) sol <- rep(NA_real_, nCol)
                mixCoef[i, ] <- sol
            } else {
                mixCoef[i, ] <- .champ_try(solve(Dmat, dvec), default = rep(NA_real_, nCol))
            }
        }
        mixCoef
    }

    ref <- .champ_cell_type_means(at)
    if (is.null(ref)) {
        .champ_stop("Cell type reference profiles are not available (ChAMPdata::",
                    if (identical(at, "27K")) "CellTypeMeans27K" else "CellTypeMeans450K",
                    "). Reinstall ChAMPdata.")
    }
    ref <- as.matrix(ref)
    pf <- .champ_probe_features(at)
    pf <- pf[rownames(beta), , drop = FALSE]
    dedup <- !duplicated(pf$Name)
    beta_dedup <- beta[dedup, , drop = FALSE]
    rownames(beta_dedup) <- pf$Name[dedup]

    common <- intersect(rownames(ref), rownames(beta_dedup))
    if (length(common) < ncol(ref) + 5L) {
        .champ_stop("Too few reference probes (", length(common), ") overlap the data.")
    }
    .champ_msg("Using ", length(common), " reference probes.")
    cellFrac <- projectWBC(beta_dedup[common, , drop = FALSE],
                           ref[common, , drop = FALSE], lessThanOne = TRUE)
    .champ_msg("Mean cell proportion estimates:")
    print(colMeans(cellFrac, na.rm = TRUE))
    if (all(is.na(cellFrac))) .champ_stop("Cell proportion estimation failed for every sample.")

    least <- which.min(colMeans(cellFrac, na.rm = TRUE))
    .champ_msg(colnames(cellFrac)[least],
               " has the smallest proportion and is used as the reference; ",
               "the other proportions are regressed out.")
    covariates <- cellFrac[, -least, drop = FALSE]
    ok <- stats::complete.cases(covariates)
    if (sum(ok) < ncol(covariates) + 2L) {
        .champ_stop("Not enough samples with complete cell proportion estimates.")
    }
    lm.o <- stats::lm(t(beta) ~ covariates)
    tmp.m <- t(stats::residuals(lm.o)) + rowMeans(beta)
    pos <- tmp.m[tmp.m > 0]
    if (any(tmp.m <= 0) && length(pos)) tmp.m[tmp.m <= 0] <- min(pos)
    below1 <- tmp.m[tmp.m < 1]
    if (any(tmp.m >= 1) && length(below1)) tmp.m[tmp.m >= 1] <- max(below1)
    dimnames(tmp.m) <- dimnames(beta)

    .champ_msg("[<<<< ChAMP.REFBASE END >>>>]")
    .champ_msg("[===========================]")
    list(CorrectedBeta = tmp.m, CellFraction = cellFrac)
}
