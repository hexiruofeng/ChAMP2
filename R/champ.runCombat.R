# =============================================================================
#  champ.runCombat -- sequential ComBat correction of several batch variables
# =============================================================================

#' Remove batch effects with ComBat
#'
#' Corrects one or more batch variables sequentially, keeping the phenotype of
#' interest (\code{variablename}) in the model so that biological variation is
#' preserved.
#'
#' @param beta beta value matrix.
#' @param pd sample sheet; must contain \code{variablename} and every entry of
#'   \code{batchname}.
#' @param variablename name of the phenotype column that must be protected.
#' @param batchname one or more batch columns to correct.
#' @param logitTrans apply ComBat on the logit scale (recommended for beta
#'   values).
#' @param verbose print progress.
#'
#' @return corrected beta value matrix.
#' @export
champ.runCombat <- function(beta = NULL, pd = NULL, variablename = "Sample_Group",
                            batchname = c("Slide"), logitTrans = TRUE,
                            verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<< CHAMP.RUNCOMBAT START >>]")
    .champ_msg("-----------------------------")
    .champ_require("sva", "champ.runCombat()")
    beta <- .champ_as_beta(beta)
    if (is.null(pd) || is.list(pd) && !is.data.frame(pd)) {
        .champ_stop("'pd' must be a data.frame with at least two factors.")
    }
    pd <- as.data.frame(pd, stringsAsFactors = FALSE)
    if (nrow(pd) != ncol(beta)) {
        .champ_stop("'pd' has ", nrow(pd), " rows but beta has ", ncol(beta), " samples.")
    }
    if (anyNA(beta)) {
        .champ_msg(sum(is.na(beta)), " NA values detected; imputing before ComBat.")
        beta <- champ_knn_impute(beta, k = 5L, strict = FALSE, verbose = FALSE)
    }
    if (!variablename %in% colnames(pd)) {
        .champ_stop("'variablename' (", variablename, ") is not a column of pd.")
    }
    batchname <- as.character(batchname)
    bad <- setdiff(batchname, colnames(pd))
    if (length(bad)) .champ_stop("These batch variables are not in pd: ", paste(bad, collapse = ", "))
    if (variablename %in% batchname) {
        .champ_stop("'variablename' and 'batchname' overlap; the phenotype cannot be corrected away.")
    }

    for (b in batchname) {
        if (length(unique(pd[[b]])) < 2L) {
            .champ_stop("Batch variable '", b, "' has a single value; nothing to correct.")
        }
        pd[[b]] <- factor(pd[[b]])
    }
    pd[[variablename]] <- if (.champ_is_numeric(pd[[variablename]])) {
        as.numeric(pd[[variablename]])
    } else {
        factor(pd[[variablename]])
    }
    if (length(unique(pd[[variablename]])) < 2L) {
        .champ_stop("'variablename' must contain at least two groups.")
    }

    beta_2 <- beta
    if (min(beta) <= 0 && isTRUE(logitTrans)) {
        .champ_msg("Values <= 0 replaced with the smallest positive value.")
        beta[beta <= 0] <- min(beta[beta > 0])
    }
    beta[beta >= 1] <- 1 - 1e-6

    innercombat <- function(beta, batch, formdf) {
        x <- if (isTRUE(logitTrans)) .champ_logit2(beta) else beta
        mod <- stats::model.matrix(formdf, data = pd)
        .champ_msg("  Model: ", paste(colnames(mod), collapse = " + "))
        combat <- .champ_try(sva::ComBat(dat = x, batch = batch, mod = mod,
                                          par.prior = TRUE), default = NULL)
        if (is.null(combat)) {
            .champ_warn("ComBat failed for this batch; the data are returned unchanged.")
            return(beta)
        }
        if (isTRUE(logitTrans)) .champ_ilogit2(combat) else combat
    }

    for (i in seq_along(batchname)) {
        .champ_msg("\n<< Start correcting ", batchname[i], " >>")
        rest <- if (i < length(batchname)) batchname[(i + 1L):length(batchname)] else character(0)
        rhs <- paste(c(variablename, rest), collapse = " + ")
        formdf <- stats::as.formula(paste("~", rhs))
        beta <- innercombat(beta, pd[[batchname[i]]], formdf)
    }

    if (is.null(beta)) {
        .champ_msg("champ.runCombat failed; returning the uncorrected data.")
        return(beta_2)
    }
    dimnames(beta) <- dimnames(beta_2)
    .champ_msg("champ.runCombat success. Corrected dataset will be returned.")
    .champ_msg("[<<< CHAMP.RUNCOMBAT END >>>]")
    .champ_msg("[===========================]")
    .champ_msg("[You may want to process champ.SVD() next.]\n")
    beta
}
