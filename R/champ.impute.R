# =============================================================================
#  champ.impute -- missing value handling
# =============================================================================

#' Impute missing beta values
#'
#' Three strategies are available:
#' \describe{
#'   \item{Combine}{remove samples and probes with too many missing values,
#'     then KNN-impute the rest (recommended for large data sets).}
#'   \item{Delete}{remove any probe with at least one missing value
#'     (recommended for small data sets).}
#'   \item{KNN}{KNN-impute every missing value.}
#' }
#'
#' @param beta beta value matrix (probes x samples).
#' @param pd sample sheet; rows are subset alongside the matrix.
#' @param method "Combine", "Delete" or "KNN".
#' @param k number of neighbours for KNN.
#' @param ProbeCutoff maximum fraction of missing values per probe (Combine).
#' @param SampleCutoff maximum fraction of missing values per sample (Combine).
#' @param verbose print progress.
#'
#' @return list(beta, pd) when pd is given, otherwise the imputed matrix.
#' @export
champ.impute <- function(beta = NULL, pd = NULL, method = "Combine", k = 5L,
                         ProbeCutoff = 0.2, SampleCutoff = 0.1, verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<< ChAMP.IMPUTE START >>>>]")
    .champ_msg("-----------------------------")
    beta <- .champ_as_beta(beta)
    method <- match.arg(method, c("Combine", "Delete", "KNN"))

    if (!anyNA(beta)) {
        .champ_msg("No missing value found; nothing to impute.")
        .champ_msg("[<<<<< ChAMP.IMPUTE END >>>>]")
        return(if (is.null(pd)) beta else list(beta = beta, pd = pd))
    }

    if (identical(method, "Combine")) {
        .champ_msg("<method>:Combine. (Suitable for Large Data Set)")
        colNA <- colSums(is.na(beta)) / nrow(beta)
        colValid <- which(colNA < SampleCutoff)
        .champ_msg("(1): ", ncol(beta) - length(colValid), " samples contain ",
                   SampleCutoff * 100, "% or more NA and will be removed.")
        if (!length(colValid)) .champ_stop("All samples would be removed; relax SampleCutoff.")
        beta <- beta[, colValid, drop = FALSE]
        if (!is.null(pd)) pd <- pd[colValid, , drop = FALSE]

        rowNA <- rowSums(is.na(beta)) / ncol(beta)
        rowValid <- which(rowNA < ProbeCutoff)
        .champ_msg("(2): ", nrow(beta) - length(rowValid), " probes contain ",
                   ProbeCutoff * 100, "% or more NA and will be removed.")
        if (!length(rowValid)) .champ_stop("All probes would be removed; relax ProbeCutoff.")
        beta <- beta[rowValid, , drop = FALSE]
        beta <- champ_knn_impute(beta, k = k, strict = FALSE, verbose = verbose)
        .champ_msg("(3): The remaining NA were imputed with KNN (k = ", k, ").")
    } else if (identical(method, "Delete")) {
        .champ_msg("<method>:Delete. (Suitable for Small Data Set)")
        colNA <- colSums(is.na(beta)) / nrow(beta)
        colValid <- which(colNA < SampleCutoff)
        .champ_msg("(1): ", ncol(beta) - length(colValid), " samples contain ",
                   SampleCutoff * 100, "% or more NA and will be removed.")
        beta <- beta[, colValid, drop = FALSE]
        if (!is.null(pd)) pd <- pd[colValid, , drop = FALSE]
        keep <- rowSums(is.na(beta)) == 0L
        .champ_msg("(2): ", sum(!keep), " probes were removed (they contained at least one NA).")
        beta <- beta[keep, , drop = FALSE]
    } else {
        .champ_msg("<method>:KNN. (Suitable for Large Data Set)")
        beta <- champ_knn_impute(beta, k = k, strict = FALSE, verbose = verbose)
        .champ_msg("All NA were imputed with KNN (k = ", k, ").")
    }

    .champ_msg("[<<<<< ChAMP.IMPUTE END >>>>]")
    .champ_msg("[===========================]")
    if (!is.null(pd)) list(beta = beta, pd = pd) else beta
}
