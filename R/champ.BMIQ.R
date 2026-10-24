# =============================================================================
#  champ.BMIQ -- Beta Mixture Quantile normalisation
#
#  Original algorithm: Teschendorff AE et al., Bioinformatics 2013.
#  This version no longer requires RPMM (archived on CRAN): the beta mixture
#  EM is provided by champ_blc(), which uses the same model and the same
#  number of EM iterations, so normalised values are reproduced to numerical
#  precision.
# =============================================================================

#' BMIQ normalisation of one sample
#'
#' @param beta.v beta values of one sample (no NA).
#' @param design.v probe design: 1 = type I, 2 = type II.
#' @param nL number of mixture components (BMIQ requires 3).
#' @param doH correct hemimethylated type II probes.
#' @param nfit number of probes used to fit the mixture (capped at the number
#'   of available probes).
#' @param th1.v initial thresholds for type I probes.
#' @param th2.v initial thresholds for type II probes (estimated when NULL).
#' @param niter number of EM iterations.
#' @param tol EM tolerance.
#' @param plots write diagnostic PDFs.
#' @param plotDir directory for the diagnostic PDFs.
#' @param sampleID sample label used in file names.
#' @param verbose print progress.
#'
#' @return list with nbeta (normalised beta), class1, class2, av1, av2, hf,
#'   th1, th2.
#' @export
champ.BMIQ <- function(beta.v, design.v, nL = 3, doH = TRUE, nfit = 10000,
                       th1.v = c(0.2, 0.75), th2.v = NULL, niter = 5,
                       tol = 0.001, plots = FALSE, plotDir = ".",
                       sampleID = 1, verbose = TRUE) {
    if (length(beta.v) != length(design.v)) {
        .champ_stop("beta.v and design.v must have the same length.")
    }
    if (!any(design.v == 1) || !any(design.v == 2)) {
        .champ_stop("design.v must contain both type I (1) and type II (2) probes.")
    }
    type1.idx <- which(design.v == 1)
    type2.idx <- which(design.v == 2)
    beta1.v <- beta.v[type1.idx]
    beta2.v <- beta.v[type2.idx]

    if (min(beta1.v) <= 0) beta1.v[beta1.v <= 0] <- min(beta1.v[beta1.v > 0])
    if (min(beta2.v) <= 0) beta2.v[beta2.v <= 0] <- min(beta2.v[beta2.v > 0])
    if (max(beta1.v) >= 1) beta1.v[beta1.v >= 1] <- max(beta1.v[beta1.v < 1])
    if (max(beta2.v) >= 1) beta2.v[beta2.v >= 1] <- max(beta2.v[beta2.v < 1])

    w0.m <- matrix(0, nrow = length(beta1.v), ncol = nL)
    w0.m[beta1.v <= th1.v[1L], 1L] <- 1
    w0.m[beta1.v > th1.v[1L] & beta1.v <= th1.v[2L], 2L] <- 1
    w0.m[beta1.v > th1.v[2L], 3L] <- 1

    if (isTRUE(verbose)) .champ_msg("Fitting EM beta mixture to type1 probes")
    set.seed(1234567)
    nfit <- min(nfit, length(beta1.v))
    rand.idx <- sample.int(length(beta1.v), nfit, replace = FALSE)
    em1.o <- champ_blc(matrix(beta1.v[rand.idx], ncol = 1), w = w0.m[rand.idx, ],
                       maxiter = niter, tol = tol, verbose = FALSE)
    subsetclass1.v <- apply(em1.o$w, 1L, which.max)
    edge_mean <- function(lo, hi, fallback) {
        if (!length(lo) || !length(hi)) return(fallback)
        out <- mean(c(max(lo), min(hi)))
        if (!is.finite(out)) fallback else out
    }
    th_1 <- edge_mean(beta1.v[rand.idx[subsetclass1.v == 1L]],
                      beta1.v[rand.idx[subsetclass1.v == 2L]], th1.v[1L])
    th_2 <- edge_mean(beta1.v[rand.idx[subsetclass1.v == 2L]],
                      beta1.v[rand.idx[subsetclass1.v == 3L]], th1.v[2L])
    if (!(th_1 < th_2)) {
        th_1 <- th1.v[1L]
        th_2 <- th1.v[2L]
    }
    subsetth1.v <- c(th_1, th_2)
    class1.v <- rep(2L, length(beta1.v))
    class1.v[beta1.v < subsetth1.v[1L]] <- 1L
    class1.v[beta1.v > subsetth1.v[2L]] <- 3L
    nth1.v <- subsetth1.v

    if (isTRUE(plots)) {
        tmpL.v <- as.vector(stats::rmultinom(1L, length(beta1.v), prob = em1.o$eta))
        tmpB.v <- unlist(lapply(seq_len(nL), function(l)
            stats::rbeta(tmpL.v[l], em1.o$a[l, 1L], em1.o$b[l, 1L])))
        grDevices::pdf(file.path(plotDir, paste0("Type1fit-", sampleID, ".pdf")),
                       width = 6, height = 4)
        graphics::plot(stats::density(beta1.v))
        d.o <- stats::density(tmpB.v)
        graphics::points(d.o$x, d.o$y, col = "green", type = "l")
        graphics::legend(x = 0.5, y = 3, legend = c("obs", "fit"),
                         fill = c("black", "green"), bty = "n")
        grDevices::dev.off()
    }

    d1U.o <- stats::density(beta1.v[class1.v == 1L])
    d1M.o <- stats::density(beta1.v[class1.v == 3L])
    mod1U <- d1U.o$x[which.max(d1U.o$y)]
    mod1M <- d1M.o$x[which.max(d1M.o$y)]
    d2U.o <- stats::density(beta2.v[beta2.v < 0.4])
    d2M.o <- stats::density(beta2.v[beta2.v > 0.6])
    mod2U <- d2U.o$x[which.max(d2U.o$y)]
    mod2M <- d2M.o$x[which.max(d2M.o$y)]

    if (is.null(th2.v)) {
        th2.v <- c(nth1.v[1L] + (mod2U - mod1U), nth1.v[2L] + (mod2M - mod1M))
    }
    w0.m <- matrix(0, nrow = length(beta2.v), ncol = nL)
    w0.m[beta2.v <= th2.v[1L], 1L] <- 1
    w0.m[beta2.v > th2.v[1L] & beta2.v <= th2.v[2L], 2L] <- 1
    w0.m[beta2.v > th2.v[2L], 3L] <- 1

    if (isTRUE(verbose)) .champ_msg("Fitting EM beta mixture to type2 probes")
    set.seed(1234567)
    nfit2 <- min(nfit, length(beta2.v))
    rand.idx <- sample.int(length(beta2.v), nfit2, replace = FALSE)
    em2.o <- champ_blc(matrix(beta2.v[rand.idx], ncol = 1), w = w0.m[rand.idx, ],
                       maxiter = niter, tol = tol, verbose = FALSE)

    if (isTRUE(plots)) {
        tmpL.v <- as.vector(stats::rmultinom(1L, length(beta2.v), prob = em2.o$eta))
        tmpB.v <- unlist(lapply(seq_len(nL), function(l)
            stats::rbeta(tmpL.v[l], em2.o$a[l, 1L], em2.o$b[l, 1L])))
        grDevices::pdf(file.path(plotDir, paste0("Type2fit-", sampleID, ".pdf")),
                       width = 6, height = 4)
        graphics::plot(stats::density(beta2.v))
        d.o <- stats::density(tmpB.v)
        graphics::points(d.o$x, d.o$y, col = "green", type = "l")
        graphics::legend(x = 0.5, y = 3, legend = c("obs", "fit"),
                         fill = c("black", "green"), bty = "n")
        grDevices::dev.off()
    }

    subsetclass2.v <- apply(em2.o$w, 1L, which.max)
    th2_1 <- edge_mean(beta2.v[rand.idx[subsetclass2.v == 1L]],
                       beta2.v[rand.idx[subsetclass2.v == 2L]], th2.v[1L])
    th2_2 <- edge_mean(beta2.v[rand.idx[subsetclass2.v == 2L]],
                       beta2.v[rand.idx[subsetclass2.v == 3L]], th2.v[2L])
    if (!(th2_1 < th2_2)) {
        th2_1 <- th2.v[1L]
        th2_2 <- th2.v[2L]
    }
    subsetth2.v <- c(th2_1, th2_2)
    class2.v <- rep(2L, length(beta2.v))
    class2.v[beta2.v < subsetth2.v[1L]] <- 1L
    class2.v[beta2.v > subsetth2.v[2L]] <- 3L

    classAV1.v <- as.numeric(em1.o$mu[, 1L])
    classAV2.v <- as.numeric(em2.o$mu[, 1L])

    if (isTRUE(verbose)) .champ_msg("Start normalising type 2 probes")
    nbeta2.v <- beta2.v

    lt <- 1L
    selU.idx <- which(class2.v == lt)
    selUR.idx <- selU.idx[beta2.v[selU.idx] > classAV2.v[lt]]
    selUL.idx <- selU.idx[beta2.v[selU.idx] < classAV2.v[lt]]
    p.v <- stats::pbeta(beta2.v[selUR.idx], em2.o$a[lt, 1L], em2.o$b[lt, 1L],
                        lower.tail = FALSE)
    nbeta2.v[selUR.idx] <- stats::qbeta(p.v, em1.o$a[lt, 1L], em1.o$b[lt, 1L],
                                        lower.tail = FALSE)
    p.v <- stats::pbeta(beta2.v[selUL.idx], em2.o$a[lt, 1L], em2.o$b[lt, 1L],
                        lower.tail = TRUE)
    nbeta2.v[selUL.idx] <- stats::qbeta(p.v, em1.o$a[lt, 1L], em1.o$b[lt, 1L],
                                        lower.tail = TRUE)

    lt <- 3L
    selM.idx <- which(class2.v == lt)
    selMR.idx <- selM.idx[beta2.v[selM.idx] > classAV2.v[lt]]
    selML.idx <- selM.idx[beta2.v[selM.idx] < classAV2.v[lt]]
    p.v <- stats::pbeta(beta2.v[selMR.idx], em2.o$a[lt, 1L], em2.o$b[lt, 1L],
                        lower.tail = FALSE)
    nbeta2.v[selMR.idx] <- stats::qbeta(p.v, em1.o$a[lt, 1L], em1.o$b[lt, 1L],
                                        lower.tail = FALSE)

    hf <- NA_real_
    if (isTRUE(doH)) {
        lt <- 2L
        selH.idx <- c(which(class2.v == lt), selML.idx)
        minH <- min(beta2.v[selH.idx])
        maxH <- max(beta2.v[selH.idx])
        deltaH <- maxH - minH
        deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
        deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])
        nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM
        nminH <- max(nbeta2.v[selU.idx]) + deltaUH
        ndeltaH <- nmaxH - nminH
        if (is.finite(deltaH) && deltaH > 0) {
            hf <- ndeltaH / deltaH
            nbeta2.v[selH.idx] <- nminH + hf * (beta2.v[selH.idx] - minH)
        }
    }

    pnbeta.v <- beta.v
    pnbeta.v[type1.idx] <- beta1.v
    pnbeta.v[type2.idx] <- nbeta2.v

    if (isTRUE(plots)) {
        d1.o <- stats::density(beta1.v)
        d2.o <- stats::density(beta2.v)
        d2n.o <- stats::density(nbeta2.v)
        ymax <- max(d2.o$y, d1.o$y, d2n.o$y)
        grDevices::pdf(file.path(plotDir, paste0("CheckBMIQ-", sampleID, ".pdf")),
                       width = 6, height = 4)
        graphics::plot(stats::density(beta2.v), type = "l", ylim = c(0, ymax),
                       xlim = c(0, 1))
        graphics::points(d1.o$x, d1.o$y, col = "red", type = "l")
        graphics::points(d2n.o$x, d2n.o$y, col = "blue", type = "l")
        graphics::legend(x = 0.5, y = ymax, legend = c("type1", "type2", "type2-BMIQ"),
                         bty = "n", fill = c("red", "black", "blue"))
        grDevices::dev.off()
    }

    if (isTRUE(verbose)) .champ_msg("Finished for sample ", sampleID)
    list(nbeta = pnbeta.v, class1 = class1.v, class2 = class2.v,
         av1 = classAV1.v, av2 = classAV2.v, hf = hf,
         th1 = nth1.v, th2 = th2.v)
}

#' Check a BMIQ normalised profile
#'
#' @param beta.v original beta values.
#' @param design.v probe design (1/2).
#' @param pnbeta.v BMIQ normalised beta values.
#' @return invisibly, a data.frame with the observed and normalised densities.
#' @export
CheckBMIQ <- function(beta.v, design.v, pnbeta.v) {
    type1.idx <- which(design.v == 1)
    type2.idx <- which(design.v == 2)
    beta1.v <- beta.v[type1.idx]
    beta2.v <- beta.v[type2.idx]
    pnbeta2.v <- pnbeta.v[type2.idx]
    d1 <- stats::density(beta1.v)
    d2 <- stats::density(beta2.v)
    d2n <- stats::density(pnbeta2.v)
    out <- data.frame(
        type = c(rep("type1", length(d1$y)), rep("type2", length(d2$y)),
                 rep("type2-BMIQ", length(d2n$y))),
        x = c(d1$x, d2$x, d2n$x),
        y = c(d1$y, d2$y, d2n$y))
    invisible(out)
}