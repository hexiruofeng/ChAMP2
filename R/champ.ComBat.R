# =============================================================================
#  champ.ComBat -- batch correction
#
#  ChAMP historically shipped its own copy of ComBat (a 2007 implementation
#  with hard-coded assumptions and a latent bug in the non-parametric branch).
#  Modern sva provides a maintained, better tested ComBat, so this wrapper
#  delegates to it while keeping the historical call signature.
# =============================================================================

#' ComBat batch correction
#'
#' @param dat data matrix (features x samples).
#' @param batch batch vector.
#' @param mod model matrix of covariates (may include a "Batch" column, which
#'   is removed because ComBat adds the batch term itself).
#' @param numCovs unused, kept for backwards compatibility.
#' @param par.prior use the parametric empirical Bayes prior.
#' @param prior.plots diagnostic plots.
#' @param mean.only only correct the mean (sva >= 3.30).
#'
#' @return the corrected matrix.
#' @export
champ.ComBat <- function(dat, batch, mod, numCovs = NULL, par.prior = TRUE,
                         prior.plots = FALSE, mean.only = FALSE) {
    .champ_require("sva", "champ.ComBat()")
    dat <- as.matrix(dat)
    if (is.data.frame(mod)) mod <- as.matrix(mod)
    if (!is.null(mod)) {
        if ("Batch" %in% colnames(mod)) {
            mod <- mod[, colnames(mod) != "Batch", drop = FALSE]
        }
        if (ncol(mod) == 0L) mod <- NULL
    }
    if (length(batch) != ncol(dat)) {
        .champ_stop("'batch' length does not match the number of columns of 'dat'.")
    }
    batch <- as.factor(batch)
    if (any(table(batch) < 2L)) {
        .champ_warn("Some batches contain a single sample; ComBat may fail.")
    }
    out <- sva::ComBat(dat = dat, batch = batch, mod = mod,
                       par.prior = par.prior, prior.plots = prior.plots,
                       mean.only = mean.only)
    dimnames(out) <- dimnames(dat)
    out
}
