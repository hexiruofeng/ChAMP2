# =============================================================================
#  champ.CNA -- copy number alterations from total intensity
# =============================================================================

#' Copy number analysis
#'
#' Quantile-normalises total intensities, computes log2 ratios against control
#' samples (ChAMP blood controls by default) and segments each sample with
#' DNAcopy.
#'
#' @param intensity total intensity matrix (probes x samples).
#' @param pheno phenotype vector.
#' @param control use control samples for the log ratio.
#' @param controlGroup name of the control group ("champCtls" uses the bundled
#'   blood controls).
#' @param sampleCNA segment every sample.
#' @param groupFreqPlots draw per-group gain/loss frequency plots.
#' @param Rplot,PDFplot plotting controls.
#' @param freqThreshold frequency threshold for the group plots.
#' @param resultsDir output directory.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param verbose print progress.
#'
#' @return list(sampleResult, groupResult).
#' @export
champ.CNA <- function(intensity = NULL, pheno = NULL, control = TRUE,
                      controlGroup = "champCtls", sampleCNA = TRUE,
                      groupFreqPlots = TRUE, Rplot = FALSE, PDFplot = TRUE,
                      freqThreshold = 0.3, resultsDir = "./CHAMP_CNA",
                      arraytype = "450K", verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<<< ChAMP.CNA START >>>>>]")
    .champ_msg("-----------------------------")
    .champ_require("DNAcopy", "champ.CNA()")
    .champ_require("preprocessCore", "champ.CNA()")
    intensity <- .champ_as_beta(intensity, name = "intensity")
    at <- .champ_arraytype(arraytype)
    pheno <- .champ_pheno(pheno, ncol(intensity))
    if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)
    .champ_msg("champ.CNA results will be saved in ", resultsDir, " .\n")

    pf <- .champ_probe_features(at)
    pf <- pf[rownames(intensity), , drop = FALSE]
    keep <- !is.na(pf$MAPINFO) & !is.na(pf$CHR) & !duplicated(pf$Name)
    intensity <- intensity[keep, , drop = FALSE]
    pf <- pf[keep, , drop = FALSE]
    if (!nrow(intensity)) .champ_stop("No annotated probe left for CNA analysis.")
    rownames(intensity) <- pf$Name

    .champ_msg("champ.CNA does not correct intensity batch effects; ",
               "use champ.runCombat() on the beta values beforehand if needed.")

    if (isTRUE(control)) {
        .champ_msg("<< Create control data >>")
        if (!identical(controlGroup, "champCtls") && !(controlGroup %in% pheno)) {
            .champ_msg("'", controlGroup, "' is not in the sample sheet; ",
                       "the ChAMP blood controls will be used instead.")
            controlGroup <- "champCtls"
        }
        if (identical(controlGroup, "champCtls")) {
            bc <- .champ_blood_controls()
            if (is.null(bc)) {
                .champ_stop("The ChAMP blood controls (ChAMPdata::champBloodCtls) ",
                            "are not available. Set control = FALSE or install ChAMPdata.")
            }
            ctlIntensity <- bc$intensity
            common_cpg <- intersect(rownames(intensity), rownames(ctlIntensity))
            if (!length(common_cpg)) .champ_stop("No common probe between the data and the blood controls.")
            intensity <- cbind(intensity[common_cpg, , drop = FALSE],
                               ctlIntensity[common_cpg, , drop = FALSE])
            pheno <- c(pheno, as.character(bc$pd$Sample_Group))
            .champ_msg("The blood control set contains ", ncol(ctlIntensity),
                       " control samples.")
        }
    }

    .champ_msg("<< Normalise intensities >>")
    intsqn <- preprocessCore::normalize.quantiles(as.matrix(intensity))
    colnames(intsqn) <- colnames(intensity)
    intsqnlog <- log2(intsqn + 1)
    intsqnlog[!is.finite(intsqnlog)] <- NA_real_

    samples <- if (isTRUE(control)) which(!pheno %in% controlGroup) else seq_along(pheno)
    if (!length(samples)) .champ_stop("No test sample left after removing the controls.")

    if (isTRUE(control)) {
        .champ_msg("<< Log2 ratio against the mean control >>")
        ref <- rowMeans(intsqnlog[, pheno %in% controlGroup, drop = FALSE], na.rm = TRUE)
    } else {
        .champ_msg("<< Log2 ratio against the mean of all samples >>")
        ref <- rowMeans(intsqnlog, na.rm = TRUE)
    }
    intsqnlogratio <- intsqnlog[, samples, drop = FALSE] - ref
    intsqnlogratio[!is.finite(intsqnlogratio)] <- 0

    ints <- data.frame(intensity, CHR = pf$CHR, MAPINFO = pf$MAPINFO,
                       stringsAsFactors = FALSE)
    ints$CHR[ints$CHR == "X"] <- "23"
    ints$CHR[ints$CHR == "Y"] <- "24"
    CHR <- suppressWarnings(as.numeric(ints$CHR))
    MAPINFO <- as.numeric(ints$MAPINFO)
    ok <- !is.na(CHR) & !is.na(MAPINFO)
    CHR <- CHR[ok]; MAPINFO <- MAPINFO[ok]
    intsqnlogratio <- intsqnlogratio[ok, , drop = FALSE]
    ints <- ints[ok, , drop = FALSE]

    sampleResult <- list()
    if (isTRUE(sampleCNA)) {
        .champ_msg("<< Segment samples >>")
        for (i in seq_len(ncol(intsqnlogratio))) {
            sid <- colnames(intsqnlogratio)[i]
            CNA.object <- DNAcopy::CNA(cbind(intsqnlogratio[, i]), CHR, MAPINFO,
                                       data.type = "logratio", sampleid = paste0(sid, "qn"))
            smoothed <- DNAcopy::smooth.CNA(CNA.object)
            seg <- DNAcopy::segment(smoothed, verbose = 0, alpha = 0.001,
                                    undo.splits = "sdundo", undo.SD = 2)
            sampleResult[[sid]] <- seg$output
            if (isTRUE(Rplot)) graphics::plot(seg, plot.type = "w", ylim = c(-6, 6))
            if (isTRUE(PDFplot)) {
                grDevices::pdf(file.path(resultsDir, paste0(sid, "qn.pdf")))
                graphics::plot(seg, plot.type = "w", ylim = c(-6, 6))
                grDevices::dev.off()
            }
        }
    }

    groupResult <- list()
    if (isTRUE(groupFreqPlots)) {
        .champ_msg("<< Segment groups >>")
        groups <- if (isTRUE(control)) setdiff(unique(pheno), controlGroup) else unique(pheno)
        for (g in groups) {
            sel <- which(pheno[samples] == g)
            if (!length(sel)) next
            data_group <- intsqnlogratio[, sel, drop = FALSE]
            group.CNA.object <- DNAcopy::CNA(data_group, CHR, MAPINFO,
                                             data.type = "logratio",
                                             sampleid = paste(colnames(data_group), g, "qn"))
            group.smoothed <- DNAcopy::smooth.CNA(group.CNA.object)
            group.seg <- DNAcopy::segment(group.smoothed, verbose = 0, alpha = 0.001,
                                          undo.splits = "sdundo", undo.SD = 2)
            groupResult[[g]] <- group.seg$output
            group.freq <- DNAcopy::glFrequency(group.seg, freqThreshold)

            innerplot <- function() {
                ord <- order(CHR, MAPINFO)
                labels_chr <- data.matrix(summary(as.factor(CHR)))
                test1 <- data.frame(labels_chr, row.names(labels_chr),
                                    stringsAsFactors = FALSE)
                test <- data.frame(chr = unique(CHR))
                colnames(test1) <- c("count", "chr")
                F1 <- merge(test, test1, by = "chr", sort = TRUE)
                if (nrow(F1) < 2L) {
                    graphics::plot.new()
                    return(invisible(NULL))
                }
                for (i in 2:nrow(F1)) F1[i, 2] <- F1[i - 1, 2] + F1[i, 2]
                F1$label <- NULL
                F1[1, 3] <- F1[1, 2] / 2
                for (i in 2:nrow(F1)) F1[i, 3] <- (F1[i, 2] + F1[i - 1, 2]) / 2
                y1 <- group.freq$gain
                y2 <- group.freq$loss
                graphics::plot(y1, type = "h", xaxt = "n", yaxt = "n", col = "green",
                               main = paste0("Frequency plot of ", g, " samples"),
                               ylim = range(-1, 1), xlab = "Chromosome number",
                               ylab = paste0("Fraction of samples with gain or loss (n=",
                                             ncol(data_group), ")"),
                               xaxs = "i", yaxs = "i")
                graphics::points(y2, type = "h", col = "red")
                x <- c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)
                graphics::axis(1, at = F1[, 2], labels = FALSE, las = 1, lty = "dotted", tck = 1)
                graphics::axis(1, at = F1[, 3], labels = F1$chr, tick = FALSE)
                graphics::axis(2, at = x, labels = x, las = 1, lty = "dotted", tck = 1)
                invisible(NULL)
            }
            if (isTRUE(Rplot)) innerplot()
            if (isTRUE(PDFplot)) {
                grDevices::pdf(file.path(resultsDir, paste0(g, "_FreqPlot.pdf")),
                               width = 10, height = 9)
                innerplot()
                grDevices::dev.off()
            }
        }
    }

    .champ_msg("[<<<<<< ChAMP.CNA END >>>>>>]")
    .champ_msg("[===========================]")
    list(sampleResult = sampleResult, groupResult = groupResult)
}
