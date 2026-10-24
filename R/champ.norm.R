# =============================================================================
#  champ.norm -- normalisation
# =============================================================================

#' Normalise methylation data
#'
#' @param beta beta value matrix.
#' @param rgSet RGChannelSet (required for SWAN and FunctionalNormalization).
#' @param mset MethylSet (required for SWAN).
#' @param resultsDir output directory for BMIQ diagnostic plots.
#' @param method "BMIQ" (default), "PBC", "SWAN" or "FunctionalNormalization".
#' @param plotBMIQ write BMIQ diagnostic PDFs.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param cores number of workers for BMIQ.
#' @param verbose print progress.
#'
#' @return normalised beta value matrix.
#' @export
champ.norm <- function(beta = NULL, rgSet = NULL, mset = NULL,
                       resultsDir = "./CHAMP_Normalization/", method = "BMIQ",
                       plotBMIQ = FALSE, arraytype = "450K", cores = 3L,
                       verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[>>>>> ChAMP.NORM START <<<<<<]")
    .champ_msg("-----------------------------")
    beta <- .champ_as_beta(beta)
    at <- .champ_arraytype(arraytype)
    method <- match.arg(method, c("BMIQ", "PBC", "SWAN", "FunctionalNormalization"))
    if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)
    .champ_msg("champ.norm results will be saved in ", resultsDir)

    if (anyNA(beta)) {
        .champ_msg(sum(is.na(beta)), " NA values found; they are imputed before normalisation.")
        beta <- champ_knn_impute(beta, k = 5L, strict = FALSE, verbose = FALSE)
    }
    if (min(beta, na.rm = TRUE) <= 0) {
        beta[beta <= 0] <- 1e-6
        .champ_msg("Values <= 0 were replaced with 1e-6.")
    }

    if (method %in% c("SWAN", "FunctionalNormalization")) {
        .champ_require("minfi", paste0("method = '", method, "'"))
        if (is.null(rgSet)) .champ_stop("'rgSet' is required for ", method, ".")
    }

    if (identical(method, "SWAN")) {
        if (is.null(mset)) .champ_stop("'mset' is required for SWAN.")
        beta.p <- minfi::getBeta(minfi::preprocessSWAN(rgSet, mset), "Illumina")
        beta.p <- beta.p[rownames(beta), , drop = FALSE]
    } else if (identical(method, "FunctionalNormalization")) {
        if (identical(at, "EPICv1")) rgSet@annotation[2L] <- "ilm10b4.hg19"
        beta.p <- minfi::getBeta(minfi::preprocessFunnorm(rgSet))[rownames(beta), , drop = FALSE]
    } else {
        pf <- .champ_probe_features(at)
        type <- pf[rownames(beta), "Type"]
        if (all(is.na(type))) .champ_stop("The annotation has no probe design type; ",
                                          method, " cannot be applied.")
        design.v <- ifelse(type == "I", 1L, 2L)
        design.v[is.na(design.v)] <- 2L
        if (!all(c(1L, 2L) %in% design.v)) {
            .champ_stop("Both type I and type II probes are needed for ", method, ".")
        }
        if (isTRUE(plotBMIQ)) dir.create(resultsDir, showWarnings = FALSE, recursive = TRUE)

        if (identical(method, "PBC")) {
            .champ_msg("<< Normalizing data with PBC Method >>")
            beta.p <- DoPBC(beta, design.v, verbose = verbose)
        } else {
            .champ_msg("<< Normalizing data with BMIQ Method >>")
            .champ_msg("BMIQ may fail for low quality samples (no beta distribution).")
            cores <- .champ_cores(cores)
            one <- function(x) {
                champ.BMIQ(beta[, x], design.v, sampleID = colnames(beta)[x],
                           plots = plotBMIQ, plotDir = resultsDir, verbose = FALSE)$nbeta
            }
            if (cores > 1L && .champ_has("parallel")) {
                .champ_msg(cores, " cores will be used for BMIQ.")
                cl <- parallel::makeCluster(cores)
                on.exit(parallel::stopCluster(cl), add = TRUE)
                res <- parallel::parLapply(cl, seq_len(ncol(beta)), one)
            } else {
                res <- lapply(seq_len(ncol(beta)), one)
            }
            beta.p <- do.call(cbind, res)
        }
    }

    beta.p <- as.matrix(beta.p)
    rownames(beta.p) <- rownames(beta)
    colnames(beta.p) <- colnames(beta)
    beta.p[beta.p < 0] <- 0
    beta.p[beta.p > 1] <- 1
    .champ_msg("[>>>>> ChAMP.NORM END <<<<<<]")
    .champ_msg("[===========================]")
    .champ_msg("[You may want to process champ.SVD() next.]\n")
    beta.p
}
