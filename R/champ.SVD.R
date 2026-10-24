# =============================================================================
#  champ.SVD -- singular value decomposition based batch effect detection
# =============================================================================

#' Random matrix theory estimate of the number of latent dimensions
#'
#' @param data.m numeric matrix (probes x samples).
#' @return list with cor, dim, estdens and thdens.
#' @export
EstDimRMT <- function(data.m) {
    data.m <- as.matrix(data.m)
    M <- data.m
    for (c in seq_len(ncol(M))) {
        s <- stats::var(data.m[, c])
        if (!is.finite(s) || s <= 0) s <- 1
        M[, c] <- (data.m[, c] - mean(data.m[, c])) / sqrt(s)
    }
    sigma2 <- stats::var(as.vector(M))
    Q <- nrow(data.m) / ncol(data.m)
    thdens <- function(Q, sigma2, ns) {
        lambdaMAX <- sigma2 * (1 + 1 / Q + 2 * sqrt(1 / Q))
        lambdaMIN <- sigma2 * (1 + 1 / Q - 2 * sqrt(1 / Q))
        delta <- lambdaMAX - lambdaMIN
        roundN <- 3
        step <- round(delta / ns, roundN)
        while (step == 0) {
            roundN <- roundN + 1
            step <- round(delta / ns, roundN)
        }
        lambda.v <- seq(lambdaMIN, lambdaMAX, by = step)
        dens.v <- (Q / (2 * pi * sigma2)) * sqrt((lambdaMAX - lambda.v) *
                                                     (lambda.v - lambdaMIN)) / lambda.v
        list(min = lambdaMIN, max = lambdaMAX, step = step,
             lambda = lambda.v, dens = dens.v)
    }
    th.o <- thdens(Q, sigma2, ncol(data.m))
    C <- 1 / nrow(M) * crossprod(M)
    eigen.o <- eigen(C, symmetric = TRUE)
    estdens.o <- stats::density(eigen.o$values, from = min(eigen.o$values),
                                to = max(eigen.o$values), cut = 0)
    intdim <- length(which(eigen.o$values > th.o$max))
    list(cor = C, dim = intdim, estdens = estdens.o, thdens = th.o)
}

#' Heatmap of SVD/phenotype associations
#' @keywords internal
.champ_svd_heatmap <- function(svdPV.m) {
    myPalette <- c("darkred", "red", "orange", "pink", "white")
    breaks.v <- c(-10000, -10, -5, -2, log10(0.05), 0)
    graphics::image(x = seq_len(nrow(svdPV.m)), y = seq_len(ncol(svdPV.m)),
                    z = log10(svdPV.m), col = myPalette, breaks = breaks.v,
                    xlab = "", ylab = "", axes = FALSE,
                    main = "Singular Value Decomposition Analysis (SVD)")
    graphics::axis(1, at = seq_len(nrow(svdPV.m)),
                   labels = paste0("PC-", seq_len(nrow(svdPV.m))), las = 2)
    suppressWarnings(graphics::axis(2, at = seq_len(ncol(svdPV.m)),
                                    labels = colnames(svdPV.m), las = 2))
    graphics::legend(x = -(nrow(svdPV.m) / 2.5), y = 3,
                     legend = c("p < 1e-10", "p < 1e-5", "p < 0.01", "p < 0.05", "p > 0.05"),
                     fill = c("darkred", "red", "orange", "pink", "white"),
                     par("usr")[2], par("usr")[4], xpd = NA)
    invisible(NULL)
}

#' SVD analysis of phenotype associations
#'
#' Computes the singular value decomposition of the centred beta matrix,
#' estimates the number of informative components with random matrix theory and
#' tests each component against every factor in \code{pd} (Kruskal-Wallis for
#' categorical variables, linear regression for numeric ones).
#'
#' @param beta beta value matrix.
#' @param rgSet optional RGChannelSet used to add control probe effects.
#' @param pd sample sheet.
#' @param RGEffect include control probe intensities as covariates.
#' @param PDFplot,Rplot produce plots.
#' @param resultsDir output directory.
#' @param verbose print progress.
#'
#' @return matrix of p-values (components x factors).
#' @export
champ.SVD <- function(beta = NULL, rgSet = NULL, pd = NULL, RGEffect = FALSE,
                      PDFplot = TRUE, Rplot = TRUE,
                      resultsDir = "./CHAMP_SVDimages/", verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<<< ChAMP.SVD START >>>>>]")
    .champ_msg("-----------------------------")
    beta <- .champ_as_beta(beta)
    if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)
    .champ_msg("champ.SVD results will be saved in ", resultsDir, " .\n")
    if (anyNA(beta)) {
        .champ_msg(sum(is.na(beta)), " NA values detected; they will be imputed for the SVD.")
        beta <- champ_knn_impute(beta, k = 5L, strict = FALSE, verbose = FALSE)
    }

    if (is.null(pd) || is.list(pd) && !is.data.frame(pd)) {
        .champ_stop("'pd' must be a data.frame or matrix with at least two factors.")
    }
    pd <- as.data.frame(pd, stringsAsFactors = FALSE)
    if (nrow(pd) != ncol(beta)) {
        .champ_stop("The pd file has ", nrow(pd), " rows but beta has ",
                    ncol(beta), " samples.")
    }
    keep_cols <- !colnames(pd) %in% c("Sample_Name", "Project", "filenames", "Basename") &
        vapply(pd, function(x) length(unique(x[!is.na(x)])) > 1L, logical(1L))
    PhenoTypes.lv <- pd[, keep_cols, drop = FALSE]
    if (!ncol(PhenoTypes.lv)) {
        .champ_stop("No factor in pd has more than one value; nothing to test.")
    }
    .champ_msg("<< Factors analysed: ", paste(colnames(PhenoTypes.lv), collapse = ", "), " >>")
    if (ncol(pd) > ncol(PhenoTypes.lv)) {
        .champ_msg("<< Factors ignored: ",
                   paste(setdiff(colnames(pd), colnames(PhenoTypes.lv)), collapse = ", "), " >>")
    }

    if (isTRUE(RGEffect)) {
        if (is.null(rgSet)) {
            .champ_msg("Control probe effect requested but no rgSet was supplied; skipping.")
        } else {
            .champ_require("minfi", "RGEffect = TRUE")
            at <- if (.champ_try(rgSet@annotation[1L], default = "") %in%
                      c("IlluminaHumanMethylationEPIC", "IlluminaHumanMethylationEPICv2"))
                "EPICv1" else "450K"
            cp <- .champ_control_probes(at)
            if (is.null(cp)) {
                .champ_msg("ControlProbes annotation unavailable; skipping the rgSet effect.")
            } else {
                dataC2.m <- as.data.frame(log2(apply(cp, 1L, function(x) {
                    if (x[3L] == "Grn") minfi::getGreen(rgSet)[x[2L], ] else minfi::getRed(rgSet)[x[2L], ]
                })))
                PhenoTypes.lv <- cbind(PhenoTypes.lv, dataC2.m)
                .champ_msg("<< Control probe intensities added to the phenotype table. >>")
            }
        }
    }

    tmp.m <- beta - rowMeans(beta, na.rm = TRUE)
    rmt.o <- EstDimRMT(tmp.m)
    svd.o <- svd(tmp.m)
    topPCA <- min(20L, max(1L, rmt.o$dim))
    svdPV.m <- matrix(NA_real_, nrow = topPCA, ncol = ncol(PhenoTypes.lv),
                      dimnames = list(paste0("PC-", seq_len(topPCA)),
                                      colnames(PhenoTypes.lv)))
    for (c in seq_len(topPCA)) {
        for (f in seq_len(ncol(PhenoTypes.lv))) {
            x <- PhenoTypes.lv[[f]]
            if (.champ_is_numeric(x)) {
                svdPV.m[c, f] <- summary(stats::lm(svd.o$v[, c] ~ x))$coeff[2, 4]
            } else {
                svdPV.m[c, f] <- stats::kruskal.test(svd.o$v[, c] ~ as.factor(x))$p.value
            }
        }
    }
    .champ_msg("<< SVD p-value matrix calculated. >>")

    if (isTRUE(Rplot)) {
        graphics::par(mar = c(5, 15, 2, 1))
        .champ_svd_heatmap(svdPV.m)
    }
    if (isTRUE(PDFplot)) {
        grDevices::pdf(file.path(resultsDir, "SVDsummary.pdf"), width = 8, height = 8)
        .champ_svd_heatmap(svdPV.m)
        grDevices::dev.off()
    }
    .champ_msg("<< SVD plot saved. >>")
    .champ_msg("[<<<<<< ChAMP.SVD END >>>>>>]")
    .champ_msg("[===========================]")
    .champ_msg("[If the batch effect is not significant, continue with champ.DMP()/champ.DMR(), ",
               "otherwise run champ.runCombat() first.]\n")
    svdPV.m
}
