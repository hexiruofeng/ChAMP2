# =============================================================================
#  champ.DMR -- differentially methylated regions
#
#  Three algorithms are available:
#    * Bumphunter  (bumphunter package)
#    * ProbeLasso  (ChAMP implementation)
#    * DMRcate     (DMRcate package when installed, otherwise the identical
#                   algorithm implemented in champ_dmrcate())
# =============================================================================

#' Probe annotation as a GRanges object
#'
#' ChAMPdata no longer ships the pre-built illumina450Gr / illuminaEPICGr
#' objects, so the ranges are built on the fly from the probe annotation.
#'
#' @param probes probe identifiers.
#' @param arraytype array type.
#' @return GRanges.
#' @keywords internal
.champ_probe_gr <- function(probes, arraytype = "450K") {
    .champ_require("GenomicRanges", "ProbeLasso DMR calling")
    pf <- .champ_probe_features_for(probes, arraytype)
    ok <- !is.na(pf$MAPINFO) & !is.na(pf$CHR)
    pf <- pf[ok, , drop = FALSE]
    gr <- GenomicRanges::GRanges(
        seqnames = pf$CHR,
        ranges = IRanges::IRanges(start = as.integer(pf$MAPINFO), width = 1L),
        name = rownames(pf),
        featureCgi = pf$feat.cgi,
        cgi = pf$cgi,
        feature = pf$feature,
        geneSymbol = pf$gene,
        ensemblID = NA_character_)
    GenomeInfoDb::sortSeqlevels(gr)
}

#' Differentially methylated regions
#'
#' @param beta beta value matrix.
#' @param pheno phenotype vector with exactly two groups (ProbeLasso and
#'   DMRcate) or two or more groups (Bumphunter).
#' @param compare.group the two groups to compare when \code{pheno} has more
#'   than two levels.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param method "Bumphunter", "ProbeLasso" or "DMRcate".
#' @param minProbes minimum number of probes per region.
#' @param adjPvalDmr adjusted p-value cutoff for reported regions.
#' @param cores number of workers.
#' @param maxGap,bumphunter parameters (Bumphunter).
#' @param cutoff,pickCutoff,smooth,smoothFunction,useWeights,permutations,B,
#'   nullMethod bumphunter tuning parameters.
#' @param meanLassoRadius,minDmrSep,minDmrSize,adjPvalProbe ProbeLasso
#'   parameters.
#' @param Rplot,PDFplot,resultsDir plotting controls.
#' @param rmSNPCH,dist,mafcut DMRcate probe filtering.
#' @param fdr,lambda,C DMRcate parameters.
#' @param dmrcate.backend "auto", "package" or "internal".
#' @param verbose print progress.
#'
#' @return named list with the DMR table (and, for ProbeLasso, the probe
#'   table).
#' @export
champ.DMR <- function(beta = NULL, pheno = NULL, compare.group = NULL,
                      arraytype = "450K", method = "Bumphunter", minProbes = 7L,
                      adjPvalDmr = 0.05, cores = 3L,
                      maxGap = 300, cutoff = NULL, pickCutoff = TRUE,
                      smooth = TRUE, smoothFunction = NULL, useWeights = FALSE,
                      permutations = NULL, B = 250L, nullMethod = "bootstrap",
                      meanLassoRadius = 375, minDmrSep = 1000, minDmrSize = 50,
                      adjPvalProbe = 0.05, Rplot = FALSE, PDFplot = FALSE,
                      resultsDir = "./CHAMP_ProbeLasso/", rmSNPCH = TRUE,
                      fdr = 0.05, dist = 2, mafcut = 0.05, lambda = 1000,
                      C = 2, dmrcate.backend = c("auto", "package", "internal"),
                      verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<<< ChAMP.DMR START >>>>>]")
    .champ_msg("-----------------------------")
    beta <- .champ_as_beta(beta)
    at <- .champ_arraytype(arraytype)
    method <- match.arg(method, c("Bumphunter", "ProbeLasso", "DMRcate"))
    backend <- match.arg(dmrcate.backend)
    pheno <- .champ_pheno(pheno, ncol(beta))
    if (anyNA(beta)) {
        .champ_msg(sum(is.na(beta)), " NA values detected; imputing before DMR calling.")
        beta <- champ_knn_impute(beta, k = 5L, strict = FALSE, verbose = FALSE)
    }
    pf <- .champ_probe_features(at)

    if (method %in% c("ProbeLasso", "DMRcate")) {
        if (length(unique(pheno)) > 2L) {
            if (is.null(compare.group) || length(compare.group) != 2L ||
                !all(compare.group %in% pheno)) {
                .champ_stop(method, " needs exactly two groups. ",
                            "Set compare.group = c('A', 'B').")
            }
            sel <- pheno %in% compare.group
            beta <- beta[, sel, drop = FALSE]
            pheno <- pheno[sel]
        }
    }

    if (identical(method, "Bumphunter")) {
        .champ_require("bumphunter", "method = 'Bumphunter'")
        .champ_msg("<< Find DMR with Bumphunter >>")
        cpg.idx <- intersect(rownames(beta), rownames(pf))
        if (length(cpg.idx) < minProbes) .champ_stop("Too few annotated probes.")
        Anno <- pf[cpg.idx, c("CHR", "MAPINFO"), drop = FALSE]
        chr_num <- suppressWarnings(as.numeric(Anno$CHR))
        chr_ord <- ifelse(is.na(chr_num), 100 + match(Anno$CHR, c("X", "Y")), chr_num)
        ord <- order(chr_ord, Anno$MAPINFO)
        Anno <- Anno[ord, , drop = FALSE]
        cpg.idx <- rownames(Anno)

        cl <- bumphunter::clusterMaker(Anno$CHR, Anno$MAPINFO, maxGap = maxGap)
        names(cl) <- cpg.idx
        keep_cluster <- names(which(table(cl) > minProbes))
        .champ_msg("Detected ", length(keep_cluster), " clusters with more than ",
                   minProbes, " probes within ", maxGap, " bp.")
        bumphunter.idx <- cpg.idx[cl %in% keep_cluster]
        if (!length(bumphunter.idx)) .champ_stop("No cluster with enough probes for Bumphunter.")

        X <- cbind(1, as.numeric(factor(pheno)) - 1)
        Beta <- beta[bumphunter.idx, , drop = FALSE]
        Beta[Beta <= 0.001] <- 0.001
        Beta[Beta >= 0.999] <- 0.999
        Y <- log2(Beta / (1 - Beta))
        if (is.null(smoothFunction)) smoothFunction <- bumphunter::loessByCluster

        Bumps <- bumphunter::bumphunter(
            Y, design = X,
            chr = Anno[bumphunter.idx, ]$CHR,
            pos = Anno[bumphunter.idx, ]$MAPINFO,
            cluster = cl[bumphunter.idx],
            cutoff = cutoff, pickCutoff = pickCutoff, smooth = smooth,
            smoothFunction = smoothFunction, useWeights = useWeights,
            permutations = permutations, verbose = verbose, B = B,
            nullMethod = nullMethod)
        DMR <- Bumps$table
        if (is.null(DMR) || !nrow(DMR)) .champ_stop("No DMR detected.")
        DMR <- DMR[DMR$p.valueArea <= adjPvalDmr, , drop = FALSE]
        .champ_msg("Bumphunter detected ", nrow(DMR), " DMRs with p <= ", adjPvalDmr, ".")
        if (!nrow(DMR)) .champ_stop("No DMR passed the p-value cutoff.")
        rownames(DMR) <- paste0("DMR_", seq_len(nrow(DMR)))
        DMR <- data.frame(DMR[, 1:3, drop = FALSE],
                          width = DMR[, 3L] - DMR[, 2L],
                          strand = "*",
                          DMR[, 4:ncol(DMR), drop = FALSE],
                          check.names = FALSE)
        colnames(DMR)[1:3] <- c("seqnames", "start", "end")
        return(list(BumphunterDMR = DMR))
    }

    if (identical(method, "ProbeLasso")) {
        .champ_require("GenomicRanges", "method = 'ProbeLasso'")
        .champ_require("S4Vectors", "method = 'ProbeLasso'")
        if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)
        .champ_msg("<< Find DMR with ProbeLasso >>")
        DMP <- champ.DMP(beta = beta, pheno = pheno, adjPVal = 1, arraytype = at,
                         verbose = verbose)
        if (length(DMP) != 1L) {
            .champ_stop("ProbeLasso needs exactly one comparison; ",
                        "run champ.DMR() separately for each pair of groups.")
        }
        DMP <- DMP[[1L]]
        if (!any(DMP$adj.P.Val < adjPvalProbe)) {
            .champ_stop("No probe is significant at adjPvalProbe = ", adjPvalProbe, ".")
        }

        gr <- .champ_probe_gr(rownames(DMP), at)
        myResultsGr <- gr[match(rownames(DMP), S4Vectors::mcols(gr)$name)]
        myResultsGr$P.Value <- DMP$P.Value[match(S4Vectors::mcols(myResultsGr)$name,
                                                 rownames(DMP))]
        myResultsGr$adj.P.Val <- DMP$adj.P.Val[match(S4Vectors::mcols(myResultsGr)$name,
                                                     rownames(DMP))]
        myResultsGr <- GenomeInfoDb::sortSeqlevels(myResultsGr)
        myResultsGr <- GenomicRanges::sort(myResultsGr, ignore.strand = TRUE)
        S4Vectors::mcols(myResultsGr)$adj.P.Val <- stats::p.adjust(
            S4Vectors::mcols(myResultsGr)$P.Value, method = "BH")

        .champ_msg("<< Get closestProbe for each probe >>")
        closestProbe <- as.data.frame(GenomicRanges::distanceToNearest(
            myResultsGr, ignore.strand = TRUE))$distance
        closestProbeSp <- split(closestProbe, S4Vectors::mcols(myResultsGr)$featureCgi)
        closestProbeSp <- closestProbeSp[vapply(closestProbeSp, length, integer(1L)) > 1L]
        if (!length(closestProbeSp)) {
            .champ_stop("ProbeLasso needs at least two probes per feature/CpG-island class.")
        }
        lassoQuantileDeviation <- abs(meanLassoRadius - rowMeans(as.data.frame(
            lapply(closestProbeSp, function(x) stats::quantile(x, (1:1000) / 1000)))))
        lassoQuantileThreshold <- which.min(lassoQuantileDeviation) / 1000
        lassoSizes <- lapply(closestProbeSp, function(x)
            stats::quantile(x, lassoQuantileThreshold, na.rm = TRUE))

        .champ_msg("<< Expand ranges for each probe >>")
        myResultsGrSp <- split(myResultsGr, S4Vectors::mcols(myResultsGr)$featureCgi)
        myResultsGrSp <- myResultsGrSp[names(myResultsGrSp) %in% names(lassoSizes)]
        lassoGr <- unlist(GenomicRanges::GRangesList(mapply(function(x, y)
            GenomicRanges::promoters(x, upstream = y, downstream = y),
            x = myResultsGrSp, y = lassoSizes[names(myResultsGrSp)],
            SIMPLIFY = FALSE)))
        myResultsSigGr <- myResultsGr[S4Vectors::mcols(myResultsGr)$adj.P.Val < adjPvalProbe]
        lassoProbeCountOverlap <- GenomicRanges::countOverlaps(lassoGr, myResultsSigGr,
                                                               ignore.strand = TRUE)
        dmrGr <- GenomicRanges::reduce(lassoGr[lassoProbeCountOverlap >= minProbes],
                                       min.gapwidth = minDmrSep, ignore.strand = TRUE)
        if (!length(dmrGr)) .champ_stop("No DMR passed the ProbeLasso lasso/minProbes filter.")
        GenomicRanges::strand(dmrGr) <- "*"
        dmrGr <- dmrGr[GenomicRanges::width(dmrGr) > minDmrSize]
        if (!length(dmrGr)) .champ_stop("No DMR larger than minDmrSize was found.")

        probeIndex <- as.data.frame(GenomicRanges::findOverlaps(dmrGr, myResultsGr))
        pValuesGr <- myResultsGr[probeIndex$subjectHits, "P.Value"]
        myBetas <- beta[match(S4Vectors::mcols(pValuesGr)$name, rownames(beta)), , drop = FALSE]
        myBetas <- split(as.data.frame(myBetas), probeIndex$queryHits)
        correl <- lapply(myBetas, function(x) stats::cor(t(x)))
        weights <- lapply(correl, function(x) 1 / apply(x^2, 1L, sum))
        dmrQP <- stats::qnorm(S4Vectors::mcols(pValuesGr)$P.Value)
        dmrQP <- split(dmrQP, probeIndex$queryHits)
        dmrQPW <- mapply("*", dmrQP, weights, SIMPLIFY = FALSE)
        dmrStat <- lapply(dmrQPW, sum)
        dmrSd <- lapply(weights, function(x) sqrt(sum(x^2)))
        dmrP <- stats::p.adjust(mapply(function(x, y) stats::pnorm(x, 0, sd = y),
                                       dmrStat, dmrSd), method = "BH")
        goodDmr <- which(dmrP < adjPvalDmr)
        if (!length(goodDmr)) .champ_stop("No DMR passed the ProbeLasso combined p-value cutoff.")
        dmrGr <- dmrGr[goodDmr]
        dmrP <- dmrP[goodDmr]
        dmrpRank <- rank(dmrP, ties.method = "min")

        probeIndex <- as.data.frame(GenomicRanges::findOverlaps(dmrGr, myResultsGr))
        dmrProbesGr <- myResultsGr[probeIndex$subjectHits]
        myBetas <- beta[match(S4Vectors::mcols(dmrProbesGr)$name, rownames(beta)), , drop = FALSE]
        myBetas <- as.data.frame(myBetas)
        dmrCoreStart <- split(GenomicRanges::start(dmrProbesGr), probeIndex$queryHits)
        dmrCoreEnd <- split(GenomicRanges::end(dmrProbesGr), probeIndex$queryHits)
        dmrCoreStart <- vapply(dmrCoreStart, min, numeric(1L))
        dmrCoreEnd <- vapply(dmrCoreEnd, max, numeric(1L))
        myBetas <- split(myBetas, probeIndex$queryHits)

        groupIndex <- pheno
        dmrGroupMeans <- do.call(rbind, lapply(myBetas, function(x)
            vapply(split(t(x), groupIndex), mean, numeric(1L))))
        colnames(dmrGroupMeans) <- paste("betaAv", colnames(dmrGroupMeans), sep = "_")
        probeGroupMeans <- lapply(myBetas, function(x)
            lapply(split(as.data.frame(t(x)), groupIndex), colMeans))
        probeGroupMeans <- do.call(rbind, lapply(probeGroupMeans, function(x)
            t(do.call(rbind, x))))
        colnames(probeGroupMeans) <- paste("betaAv", colnames(probeGroupMeans), sep = "_")

        myDmrProbesGr <- myResultsGr[probeIndex$subjectHits]
        myDmrProbesGr <- GenomicRanges::GRanges(
            seqnames = GenomicRanges::seqnames(myDmrProbesGr),
            ranges = IRanges::IRanges(start = GenomicRanges::start(myDmrProbesGr),
                                      end = GenomicRanges::end(myDmrProbesGr)))
        S4Vectors::mcols(myDmrProbesGr) <- S4Vectors::DataFrame(
            cbind(as.data.frame(S4Vectors::mcols(myDmrProbesGr)), probeGroupMeans))
        myDmrProbesGr$dmrNo <- probeIndex$queryHits
        myDmrProbesGr$dmrP <- dmrP[probeIndex$queryHits]
        myDmrProbesGr$dmrpRank <- dmrpRank[probeIndex$queryHits]
        myDmrProbesGr$dmrChrom <- GenomicRanges::seqnames(dmrGr[probeIndex$queryHits])
        myDmrProbesGr$dmrStart <- GenomicRanges::start(dmrGr[probeIndex$queryHits])
        myDmrProbesGr$dmrEnd <- GenomicRanges::end(dmrGr[probeIndex$queryHits])
        myDmrProbesGr$dmrSize <- GenomicRanges::width(dmrGr[probeIndex$queryHits])
        myDmrProbesGr$dmrCoreStart <- dmrCoreStart[probeIndex$queryHits]
        myDmrProbesGr$dmrCoreEnd <- dmrCoreEnd[probeIndex$queryHits]
        myDmrProbesGr$dmrCoreSize <- myDmrProbesGr$dmrCoreEnd - myDmrProbesGr$dmrCoreStart + 1

        myDmrGr <- dmrGr
        myDmrGr$dmrNo <- unique(probeIndex$queryHits)
        myDmrGr$dmrP <- dmrP
        myDmrGr$dmrpRank <- dmrpRank
        myDmrGr$dmrChrom <- GenomicRanges::seqnames(dmrGr)
        myDmrGr$dmrStart <- GenomicRanges::start(dmrGr)
        myDmrGr$dmrEnd <- GenomicRanges::end(dmrGr)
        myDmrGr$dmrSize <- GenomicRanges::width(dmrGr)
        myDmrGr$dmrCoreStart <- dmrCoreStart
        myDmrGr$dmrCoreEnd <- dmrCoreEnd
        myDmrGr$dmrCoreSize <- myDmrGr$dmrCoreEnd - myDmrGr$dmrCoreStart + 1
        genes <- split(as.data.frame(S4Vectors::mcols(myResultsGr))[probeIndex$subjectHits,
                                                                     c("geneSymbol", "name")],
                       probeIndex$queryHits)
        myDmrGr$geneSymbol <- vapply(genes, function(x)
            paste(unique(unlist(strsplit(x$geneSymbol, ";"))), collapse = ";"), character(1L))
        myDmrGr <- GenomicRanges::GRanges(
            seqnames = GenomicRanges::seqnames(myDmrGr),
            ranges = IRanges::ranges(myDmrGr))
        S4Vectors::mcols(myDmrGr) <- S4Vectors::DataFrame(
            cbind(as.data.frame(S4Vectors::mcols(myDmrGr)), dmrGroupMeans))

        if (isTRUE(Rplot) || isTRUE(PDFplot)) {
            interplot <- function(lassoSizes) {
                cgiNames <- levels(factor(S4Vectors::mcols(myResultsGr)$cgi))
                featureNames <- levels(factor(S4Vectors::mcols(myResultsGr)$feature))
                if (length(cgiNames) < 4L || length(featureNames) < 7L) {
                    graphics::plot.new()
                    graphics::title("lasso radius plot unavailable for this annotation")
                    return(invisible(NULL))
                }
                sfo <- c(6, 7, 3, 1, 4, 2, 5)
                scgio <- c(1, 4, 3, 2)
                sfcgio <- rep((sfo - 1) * 4, each = 4) + rep(scgio, 7)
                lassoSizes <- round(unlist(lassoSizes))
                graphics::par(mar = c(7, 4, 4, 3) + 0.5)
                graphics::plot(c(1, 28),
                               y = c(range(0.3 * sqrt(lassoSizes))[1] * 0.8,
                                     range(0.3 * sqrt(lassoSizes))[2] * 1.2),
                               type = "n", xaxt = "n", xlab = "", yaxt = "n",
                               ylab = "lasso radius [bp]",
                               main = paste0("lasso quantile = ", round(lassoQuantileThreshold, 2),
                                             "\nmean lasso radius = ", meanLassoRadius, "bp"),
                               bty = "n")
                graphics::segments(1:28, rep(0, 28), 1:28,
                                   0.3 * sqrt(lassoSizes[sfcgio]), lty = 3, col = "grey")
                graphics::points(1:28, 0.3 * sqrt(lassoSizes[sfcgio]),
                                 pch = 16, cex = 0.3 * sqrt(lassoSizes[sfcgio]),
                                 col = rep(grDevices::rainbow(7, alpha = 0.5)[sfo], each = 4))
                graphics::text(1:28, 0.3 * sqrt(lassoSizes[sfcgio]), lassoSizes[sfcgio],
                               pos = 3, cex = 0.8)
                graphics::axis(1, at = 1:28, labels = rep(cgiNames[scgio], 7), las = 2)
                graphics::par(xpd = TRUE)
                graphics::segments(seq(1, 28, 4), rep(-2.5, 7), seq(4, 28, 4), rep(-2.5, 7))
                graphics::mtext(text = featureNames[sfo], side = 1, at = seq(2.5, 28, 4),
                                line = 5.5, las = 1, cex.axis = 1)
                graphics::axis(2, at = c(0, max(0.3 * sqrt(lassoSizes))), labels = FALSE)
                invisible(NULL)
            }
            if (isTRUE(Rplot)) interplot(lassoSizes)
            if (isTRUE(PDFplot)) {
                grDevices::pdf(file.path(resultsDir, "myLassos.pdf"), width = 9, height = 9)
                interplot(lassoSizes)
                grDevices::dev.off()
            }
        }

        DMRProbes <- as.data.frame(myDmrProbesGr)
        keep_cols <- which(colnames(DMRProbes) == "P.Value"):which(colnames(DMRProbes) == "dmrNo")
        DMRProbes <- data.frame(pf[rownames(DMRProbes), , drop = FALSE],
                                DMRProbes[, keep_cols, drop = FALSE], check.names = FALSE)
        DMRProbes <- split(DMRProbes, DMRProbes$dmrNo)
        DMR <- as.data.frame(myDmrGr)
        if (identical(at, "EPICv1") || identical(at, "EPICv2")) {
            DMR[, 1L] <- paste0("chr", DMR[, 1L])
        }
        rownames(DMR) <- paste0("DMR_", DMR$dmrNo)
        names(DMRProbes) <- rownames(DMR)
        .champ_msg("ProbeLasso detected ", nrow(DMR), " DMRs.")
        return(list(ProbeLassoDMR = DMR, ProbeLassoDMRProbes = DMRProbes))
    }

    ## ---- DMRcate ---------------------------------------------------------
    .champ_msg("<< Find DMR with DMRcate >>")
    use_pkg <- identical(backend, "package") ||
        (identical(backend, "auto") && .champ_has("DMRcate") &&
             !identical(at, "EPICv2") && .champ_has("minfi"))
    DMR <- NULL
    if (use_pkg) {
        .champ_require("DMRcate", "method = 'DMRcate'")
        DMR <- .champ_try({
            myMs <- .champ_logit2(beta)
            if (isTRUE(rmSNPCH)) {
                myMs <- DMRcate::rmSNPandCH(myMs, dist = dist, mafcut = mafcut)
            }
            design <- stats::model.matrix(~factor(pheno))
            ann <- DMRcate::cpg.annotate(
                datatype = "array", object = myMs, what = "M",
                arraytype = if (identical(at, "450K")) "450K" else "EPICv1",
                design = design, coef = ncol(design),
                analysis.type = "differential", fdr = fdr)
            res <- DMRcate::dmrcate(ann, min.cpgs = minProbes, lambda = lambda, C = C)
            .champ_dmr_from_coords(res, pf)
        }, default = NULL)
        if (is.null(DMR)) {
            .champ_msg("The DMRcate package backend failed; using the built-in implementation.")
        }
    }
    if (is.null(DMR)) {
        DMR <- champ_dmrcate(beta = beta, pheno = pheno, arraytype = at,
                             fdr = fdr, min.cpgs = minProbes, lambda = lambda,
                             C = C, verbose = verbose)
    }
    if (is.null(DMR) || !nrow(DMR)) .champ_stop("No DMR detected.")
    .champ_msg("DMRcate detected ", nrow(DMR), " DMRs.")
    if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)
    list(DMRcateDMR = DMR)
}

#' Convert a DMRcate DMResults object into the ChAMP DMR table
#' @keywords internal
.champ_dmr_from_coords <- function(res, pf) {
    coord <- res@coord
    chr <- sub(":.*$", "", coord)
    rng <- sub("^.*:", "", coord)
    start <- as.integer(sub("-.*$", "", rng))
    end <- as.integer(sub("^.*-", "", rng))
    out <- data.frame(
        seqnames = chr,
        start = start,
        end = end,
        width = end - start + 1L,
        strand = "*",
        no.cpgs = res@no.cpgs,
        min_smoothed_fdr = res@min_smoothed_fdr,
        Stouffer = res@Stouffer,
        HMFDR = res@HMFDR,
        Fisher = res@Fisher,
        maxdiff = res@maxdiff,
        meandiff = res@meandiff,
        stringsAsFactors = FALSE)
    out$overlapping.genes <- .champ_genes_in_regions(out, pf)
    rownames(out) <- paste0("DMR_", seq_len(nrow(out)))
    out
}