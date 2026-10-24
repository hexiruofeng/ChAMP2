# =============================================================================
#  champ.QC -- quality control plots, implemented without minfi
# =============================================================================

#' MDS plot of the most variable probes
#' @keywords internal
.champ_mds_plot <- function(beta, pheno, numPositions = 1000, main = "MDS plot",
                            legendPos = "topright") {
    v <- .champ_row_vars(beta)
    n <- min(numPositions, length(v))
    idx <- order(v, decreasing = TRUE)[seq_len(n)]
    d <- stats::dist(t(beta[idx, , drop = FALSE]))
    fit <- stats::cmdscale(d, k = 2)
    groups <- factor(pheno)
    cols <- grDevices::rainbow(length(levels(groups)))[as.integer(groups)]
    graphics::plot(fit[, 1L], fit[, 2L], pch = 19, col = cols,
                   xlab = "Dimension 1", ylab = "Dimension 2", main = main)
    graphics::legend(legendPos, legend = levels(groups),
                     col = grDevices::rainbow(length(levels(groups))), pch = 19, bty = "n")
    invisible(fit)
}

#' Density plot of beta values
#' @keywords internal
.champ_density_plot <- function(beta, pheno, main = "Density plot",
                                xlab = "Beta value") {
    groups <- factor(pheno)
    cols <- grDevices::rainbow(length(levels(groups)))
    d <- stats::density(beta[, 1L], na.rm = TRUE)
    graphics::plot(d, type = "n", main = main, xlab = xlab, ylab = "Density",
                   ylim = c(0, max(d$y) * 1.1))
    for (j in seq_len(ncol(beta))) {
        dj <- stats::density(beta[, j], na.rm = TRUE)
        graphics::lines(dj, col = cols[as.integer(groups[j])])
    }
    graphics::legend("topright", legend = levels(groups), col = cols, lwd = 1, bty = "n")
    invisible(NULL)
}

#' Sample dendrogram
#' @keywords internal
.champ_dendrogram <- function(beta, pheno, main = "Sample clustering") {
    hc <- stats::hclust(stats::dist(t(beta)))
    groups <- factor(pheno)
    cols <- grDevices::rainbow(length(levels(groups)))
    if (.champ_has("dendextend")) {
        dend <- stats::as.dendrogram(hc)
        dend <- dendextend::set(dend, "labels_cex", 0.8)
        dendextend::labels_colors(dend) <- cols[as.integer(groups[stats::order.dendrogram(dend)])]
        graphics::plot(dend, center = TRUE, main = main)
    } else {
        graphics::plot(hc, main = main, xlab = "", sub = "")
    }
    graphics::legend("topright", legend = levels(groups), fill = cols, bty = "n")
    invisible(hc)
}

#' Quality control plots
#'
#' Produces MDS, density and dendrogram plots before and/or after
#' normalisation.  The plotting helpers are self-contained (no minfi needed).
#'
#' @param beta beta value matrix.
#' @param pheno phenotype vector.
#' @param mdsPlot,densityPlot,dendrogram which plots to produce.
#' @param PDFplot write PDF files into resultsDir.
#' @param Rplot draw on the current device.
#' @param Feature.sel "None" or "SVD" probe selection for the dendrogram.
#' @param resultsDir output directory.
#' @param verbose print progress.
#'
#' @return invisibly, the hclust object when a dendrogram is produced.
#' @export
champ.QC <- function(beta = NULL, pheno = NULL, mdsPlot = TRUE, densityPlot = TRUE,
                     dendrogram = TRUE, PDFplot = TRUE, Rplot = TRUE,
                     Feature.sel = "None", resultsDir = "./CHAMP_QCimages/",
                     verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<<< ChAMP.QC START >>>>>>]")
    .champ_msg("-----------------------------")
    beta <- .champ_as_beta(beta)
    if (is.null(pheno)) pheno <- colnames(beta)
    pheno <- .champ_pheno(pheno, ncol(beta))
    if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)
    .champ_msg("champ.QC results will be saved in ", resultsDir)
    .champ_msg("[QC plots will be produced with ", nrow(beta), " probes and ",
               ncol(beta), " samples.]\n")

    if (anyNA(beta)) .champ_msg(length(which(is.na(beta))),
                                " NA values are present; plots use complete cases.")
    if (min(beta, na.rm = TRUE) <= 0) {
        beta[beta <= 0] <- min(beta[beta > 0], na.rm = TRUE)
        .champ_msg("Values <= 0 were replaced with the smallest positive value.")
    }
    .champ_msg("<< Prepare Data Over. >>")

    hc <- NULL
    if (isTRUE(mdsPlot)) {
        if (isTRUE(Rplot)) {
            .champ_mds_plot(beta, pheno)
        }
        if (isTRUE(PDFplot)) {
            grDevices::pdf(file.path(resultsDir, "raw_mdsPlot.pdf"), width = 6, height = 4)
            .champ_mds_plot(beta, pheno)
            grDevices::dev.off()
        }
        .champ_msg("<< MDS plot done. >>\n")
    }

    if (isTRUE(densityPlot)) {
        main <- paste0("Density plot of raw data (", nrow(beta), " probes)")
        if (isTRUE(Rplot)) .champ_density_plot(beta, pheno, main = main)
        if (isTRUE(PDFplot)) {
            grDevices::pdf(file.path(resultsDir, "raw_densityPlot.pdf"), width = 6, height = 4)
            .champ_density_plot(beta, pheno, main = main)
            grDevices::dev.off()
        }
        .champ_msg("<< Density plot done. >>\n")
    }

    if (isTRUE(dendrogram)) {
        b <- beta
        Feature.sel <- match.arg(Feature.sel, c("None", "SVD"))
        if (identical(Feature.sel, "SVD")) {
            .champ_msg("<Dendrogram Feature Selection Method>: top SVD components.")
            s <- svd(b - rowMeans(b, na.rm = TRUE))
            ncomp <- max(2L, min(10L, sum(s$d > 1e-8)))
            b <- t(s$v[, seq_len(ncomp), drop = FALSE])
            colnames(b) <- colnames(beta)
        } else {
            .champ_msg("<Dendrogram Feature Selection Method>: no selection.")
        }
        main <- paste0("All samples before normalization (", nrow(beta), " probes)")
        if (isTRUE(Rplot)) hc <- .champ_dendrogram(b, pheno, main = main)
        if (isTRUE(PDFplot)) {
            grDevices::pdf(file.path(resultsDir, "raw_SampleCluster.pdf"),
                           width = max(6, floor(log(ncol(beta)) * 3)),
                           height = max(5, floor(log(ncol(beta)) * 2)))
            hc <- .champ_dendrogram(b, pheno, main = main)
            grDevices::dev.off()
        }
        .champ_msg("<< Dendrogram done. >>\n")
    }

    .champ_msg("[<<<<<< ChAMP.QC END >>>>>>>]")
    .champ_msg("[===========================]")
    .champ_msg("[You may want to process champ.norm() next.]\n")
    invisible(hc)
}