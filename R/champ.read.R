# =============================================================================
#  champ.read / champ.colname -- small helpers kept for compatibility
# =============================================================================

#' Read a beta matrix and a sample sheet from delimited files
#'
#' @param betaFile tab separated beta matrix file.
#' @param sampleSheet tab separated sample sheet.
#' @param resultsDir output directory (created when missing).
#' @param sep field separator.
#' @param verbose print progress.
#'
#' @return list(beta, pd).
#' @export
champ.read <- function(betaFile = "beta.txt", sampleSheet = "sampleSheet.txt",
                       resultsDir = NULL, sep = "\t", verbose = TRUE) {
    beta <- utils::read.table(betaFile, header = TRUE, sep = sep,
                              check.names = FALSE)
    pd <- utils::read.table(sampleSheet, header = TRUE, sep = sep,
                            stringsAsFactors = FALSE)
    if (!is.null(resultsDir)) {
        if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)
    }
    list(beta = beta, pd = pd)
}

#' Rename the columns of a beta matrix with Sample_Name
#'
#' @param beta beta value matrix.
#' @param pd sample sheet containing Sample_Name.
#' @return the beta matrix with renamed columns.
#' @export
champ.colname <- function(beta, pd) {
    if (!"Sample_Name" %in% colnames(pd)) {
        .champ_stop("'pd' must contain a Sample_Name column.")
    }
    if (length(pd$Sample_Name) < ncol(beta)) {
        .champ_stop("'pd' has fewer Sample_Name entries than beta columns.")
    }
    colnames(beta) <- as.character(pd$Sample_Name[seq_len(ncol(beta))])
    beta
}
