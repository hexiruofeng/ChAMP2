# =============================================================================
#  champ.filter -- probe / sample filtering for beta (or M) matrices
# =============================================================================

#' Filter a methylation matrix
#'
#' Applies ChAMP's standard filtering steps to a beta (or M) value matrix:
#' detection p-value filtering, optional KNN imputation, bead count filtering,
#' non-CpG removal, SNP masking, multi-hit removal and sex chromosome removal.
#'
#' @param beta beta value matrix (probes x samples); can be NULL when M is
#'   given.
#' @param M M value matrix; can be NULL when beta is given.
#' @param pd sample sheet.
#' @param intensity,Meth,UnMeth optional matrices that are filtered alongside.
#' @param detP detection p-value matrix.
#' @param beadcount bead count matrix.
#' @param autoimpute impute remaining missing values (KNN, k = 5).
#' @param filterDetP,filterBeads,filterNoCG,filterSNPs,filterMultiHit,filterXY
#'   toggles for the individual filters.
#' @param ProbeCutoff fraction of samples that may fail for a probe (0 removes
#'   any probe failing in one sample).
#' @param SampleCutoff fraction of failed probes above which a sample is
#'   removed.
#' @param detPcut detection p-value threshold.
#' @param beadCutoff fraction of samples with bead count < 3 above which a
#'   probe is removed.
#' @param population SNP mask population (e.g. "EUR", "AFR"); NULL uses the
#'   general mask.
#' @param fixOutlier replace exact 0/1 beta values.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param verbose print progress.
#'
#' @return list of the filtered matrices plus pd.
#' @export
champ.filter <- function(beta = NULL, M = NULL, pd = NULL, intensity = NULL,
                         Meth = NULL, UnMeth = NULL, detP = NULL,
                         beadcount = NULL, autoimpute = TRUE,
                         filterDetP = TRUE, ProbeCutoff = 0, SampleCutoff = 0.1,
                         detPcut = 0.01, filterBeads = TRUE, beadCutoff = 0.05,
                         filterNoCG = TRUE, filterSNPs = TRUE, population = NULL,
                         filterMultiHit = TRUE, filterXY = TRUE, fixOutlier = TRUE,
                         arraytype = "450K", verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<< ChAMP.FILTER START >>>>>]")
    .champ_msg("-----------------------------")
    at <- .champ_arraytype(arraytype)

    .champ_msg("\n[ Section 1:  Check Input Start ]")
    Objects <- list(beta = beta, M = M, intensity = intensity,
                    Meth = Meth, UnMeth = UnMeth)
    Objects <- Objects[!vapply(Objects, is.null, logical(1L))]
    if (!length(Objects)) .champ_stop("  At least one data matrix is required.")
    .champ_msg("  You have inputed ", paste(names(Objects), collapse = ","), " for Analysis.")
    if (length(unique(lapply(Objects, rownames))) != 1L ||
        length(unique(lapply(Objects, colnames))) != 1L) {
        .champ_stop("  Your objects have different rownames/colnames. ",
                    "Please make sure they are in accord.")
    }
    Objects <- lapply(Objects, function(x) {
        x <- as.matrix(x)
        storage.mode(x) <- "double"
        x
    })

    FilterOption <- list(filterDetP = isTRUE(filterDetP),
                         autoimpute = isTRUE(autoimpute),
                         filterBeads = isTRUE(filterBeads),
                         filterNoCG = isTRUE(filterNoCG),
                         filterSNPs = isTRUE(filterSNPs),
                         filterMultiHit = isTRUE(filterMultiHit),
                         filterXY = isTRUE(filterXY))

    if (!is.null(pd)) {
        .champ_msg("\n  pd file provided, checking if it is in accord with the data...")
        if (nrow(pd) != ncol(Objects[[1L]])) {
            .champ_stop("    pd file and data matrix have different dimensions.")
        }
        if ("Sample_Name" %in% names(pd) &&
            !identical(as.character(pd$Sample_Name), colnames(Objects[[1L]]))) {
            .champ_stop("    Your pd file's Sample_Name differs from the matrix colnames.")
        }
        .champ_msg("    pd file check success.")
    }

    if (isTRUE(FilterOption$filterDetP)) {
        if (is.null(detP)) {
            .champ_msg("    !!! detP not supplied; filterDetP is set to FALSE.")
            FilterOption$filterDetP <- FALSE
        } else if (!identical(dim(detP), dim(Objects[[1L]]))) {
            .champ_msg("    !!! detP does not match the data matrix; filterDetP is set to FALSE.")
            FilterOption$filterDetP <- FALSE
        } else {
            detP <- detP[rownames(Objects[[1L]]), colnames(Objects[[1L]]), drop = FALSE]
        }
    }
    if (isTRUE(FilterOption$filterBeads)) {
        if (is.null(beadcount)) {
            .champ_msg("    !!! beadcount not supplied; filterBeads is set to FALSE.")
            FilterOption$filterBeads <- FALSE
        } else if (!identical(dim(beadcount), dim(Objects[[1L]]))) {
            .champ_msg("    !!! beadcount does not match the data matrix; filterBeads is set to FALSE.")
            FilterOption$filterBeads <- FALSE
        } else {
            beadcount <- beadcount[rownames(Objects[[1L]]), colnames(Objects[[1L]]), drop = FALSE]
        }
    }
    if (isTRUE(FilterOption$autoimpute)) {
        if (!any(c("beta", "M") %in% names(Objects))) {
            .champ_msg("    !!! beta or M matrix is required for imputation; autoimpute is set to FALSE.")
            FilterOption$autoimpute <- FALSE
        } else if (!(ProbeCutoff > 0 && !is.null(detP))) {
            .champ_msg("    !!! ProbeCutoff is 0, so there is nothing to impute; autoimpute is set to FALSE.")
            FilterOption$autoimpute <- FALSE
        }
    }
    if (sum(unlist(FilterOption)) == 0L) {
        .champ_stop("  No filtering would be done. Please check your input parameters.")
    }
    .champ_msg("  Checking Finished: ",
               paste(names(which(unlist(FilterOption))), collapse = ","), " will be done.")
    .champ_msg("[ Section 1: Check Input Done ]")

    .champ_msg("\n\n[ Section 2: Filtering Start >>")
    RemainSample <- seq_len(ncol(Objects[[1L]]))

    if (FilterOption$filterDetP) {
        .champ_msg("\n  Filtering Detect P value Start")
        numfail <- matrix(colMeans(detP >= detPcut))
        rownames(numfail) <- colnames(detP)
        colnames(numfail) <- "Failed CpG Fraction."
        print(numfail)
        RemainSample <- which(numfail < SampleCutoff)
        if (any(numfail >= SampleCutoff)) {
            .champ_msg("\n    The detSamplecut parameter is : ", SampleCutoff)
            .champ_msg("    Samples : ",
                       paste(rownames(numfail)[numfail >= SampleCutoff], collapse = ","),
                       " will be deleted.")
            .champ_msg("    There are ", length(RemainSample), " samples remained for analysis.")
            if (!length(RemainSample)) .champ_stop("    All samples were removed by detection p-value filtering.")
            Objects <- lapply(Objects, function(x) x[, RemainSample, drop = FALSE])
            detP <- detP[, RemainSample, drop = FALSE]
            if (!is.null(beadcount)) beadcount <- beadcount[, RemainSample, drop = FALSE]
        }

        RemainProbe <- rowSums(detP > detPcut) <= ProbeCutoff * length(RemainSample)
        .champ_msg("\n    Filtering probes with a detection p-value above ", detPcut,
                   " in ", if (ProbeCutoff == 0) "one or more samples" else
                       paste0("at least ", ProbeCutoff * 100, "% of samples"), ".")
        .champ_msg("    Removing ", sum(!RemainProbe), " probes.")
        Objects <- lapply(Objects, function(x) x[RemainProbe, , drop = FALSE])
        detP <- detP[RemainProbe, , drop = FALSE]
        if (!is.null(beadcount)) beadcount <- beadcount[RemainProbe, , drop = FALSE]
        if (sum(detP > detPcut) > 0) {
            .champ_msg("    There are still ", sum(detP > detPcut),
                       " failed probes in your data set; imputation is recommended.")
        }
    }

    if (FilterOption$autoimpute) {
        .champ_msg("\n  Autoimpute Start")
        nfail <- sum(detP > detPcut)
        if (nfail == 0L) {
            .champ_msg("    No failed probes remain; no imputation is needed.")
        } else {
            .champ_msg("    There are ", nfail, " failed probes. KNN imputation (k = 5) starts.")
            for (nm in intersect(c("beta", "M"), names(Objects))) {
                .champ_msg("    Doing imputation on the ", nm, " matrix.")
                Objects[[nm]][detP > detPcut] <- NA_real_
                Objects[[nm]] <- champ_knn_impute(Objects[[nm]], k = 5L,
                                                  strict = FALSE, verbose = verbose)
            }
        }
    }

    if (FilterOption$filterBeads) {
        .champ_msg("\n  Filtering BeadCount Start")
        RemainProbe <- rowSums(is.na(beadcount)) < beadCutoff * ncol(beadcount)
        .champ_msg("    Filtering probes with a beadcount < 3 in at least ",
                   beadCutoff * 100, "% of samples.")
        .champ_msg("    Removing ", sum(!RemainProbe), " probes.")
        Objects <- lapply(Objects, function(x) x[RemainProbe, , drop = FALSE])
        detP <- detP[RemainProbe, , drop = FALSE]
        beadcount <- beadcount[RemainProbe, , drop = FALSE]
    }

    if (FilterOption$filterNoCG) {
        .champ_msg("\n  Filtering NoCG Start")
        RemainProbe <- substr(rownames(Objects[[1L]]), 1L, 2L) == "cg"
        .champ_msg("    Only keeping CpGs, removing ", sum(!RemainProbe), " probes.")
        Objects <- lapply(Objects, function(x) x[RemainProbe, , drop = FALSE])
        detP <- detP[RemainProbe, , drop = FALSE]
        if (!is.null(beadcount)) beadcount <- beadcount[RemainProbe, , drop = FALSE]
    }

    if (FilterOption$filterSNPs) {
        .champ_msg("\n  Filtering SNPs Start")
        maskname <- .champ_mask_probes(at, population)
        if (!length(maskname)) {
            .champ_msg("    No SNP mask is available for this array; SNP filtering skipped.")
        } else {
            RemainProbe <- !rownames(Objects[[1L]]) %in% maskname
            .champ_msg("    Removing ", sum(!RemainProbe), " probes.")
            Objects <- lapply(Objects, function(x) x[RemainProbe, , drop = FALSE])
            detP <- detP[RemainProbe, , drop = FALSE]
            if (!is.null(beadcount)) beadcount <- beadcount[RemainProbe, , drop = FALSE]
        }
    }

    if (FilterOption$filterMultiHit) {
        .champ_msg("\n  Filtering MultiHit Start")
        mh <- .champ_multi_hit()
        if (is.null(mh)) {
            .champ_msg("    The multi-hit probe list is not available; filtering skipped.")
        } else {
            RemainProbe <- !rownames(Objects[[1L]]) %in% mh
            .champ_msg("    Removing ", sum(!RemainProbe), " probes.")
            Objects <- lapply(Objects, function(x) x[RemainProbe, , drop = FALSE])
            detP <- detP[RemainProbe, , drop = FALSE]
            if (!is.null(beadcount)) beadcount <- beadcount[RemainProbe, , drop = FALSE]
        }
    }

    if (FilterOption$filterXY) {
        .champ_msg("\n  Filtering XY Start")
        pf <- .champ_probe_features(at)
        RemainProbe <- rownames(Objects[[1L]]) %in% rownames(pf)[!pf$CHR %in% c("X", "Y")]
        .champ_msg("    Filtering probes on X/Y, removing ", sum(!RemainProbe), " probes.")
        Objects <- lapply(Objects, function(x) x[RemainProbe, , drop = FALSE])
        detP <- detP[RemainProbe, , drop = FALSE]
        if (!is.null(beadcount)) beadcount <- beadcount[RemainProbe, , drop = FALSE]
    }

    if (!is.null(pd)) {
        .champ_msg("\n  Updating PD file")
        if (!FilterOption$filterDetP) {
            RemainSample <- seq_len(nrow(pd))
        }
        pd <- pd[RemainSample, , drop = FALSE]
        rownames(pd) <- NULL
        Objects <- c(Objects, list(pd = pd))
    }

    if (isTRUE(fixOutlier) && "beta" %in% names(Objects)) {
        .champ_msg("\n  Fixing Outliers Start")
        v <- Objects$beta
        pos <- v[v > 0]
        if (any(v <= 0) && length(pos)) {
            .champ_msg("    Replacing values <= 0 with the smallest positive value.")
            v[v <= 0] <- min(pos)
        }
        below <- v[v < 1]
        if (any(v >= 1) && length(below)) {
            .champ_msg("    Replacing values >= 1 with the largest value below 1.")
            v[v >= 1] <- max(below)
        }
        Objects$beta <- v
    }
    .champ_msg("[ Section 2: Filtering Done ]")

    if (!length(Objects)) .champ_stop("  All probes were removed by filtering.")
    .champ_msg("\n All filterings are Done, now you have ", nrow(Objects[[1L]]),
               " probes and ", ncol(Objects[[1L]]), " samples.\n")
    .champ_msg("[<<<<< ChAMP.FILTER END >>>>>>]")
    .champ_msg("[===========================]")
    .champ_msg("[You may want to process champ.QC() next.]\n")
    attr(Objects, "arraytype") <- at
    Objects
}
