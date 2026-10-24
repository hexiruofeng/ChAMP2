# =============================================================================
#  champ.Block -- differentially methylated blocks (large-scale regions)
# =============================================================================

#' Differentially methylated blocks
#'
#' Probes are aggregated into CpG-island / open-sea / shore-shelf clusters and
#' blocks of spatially adjacent clusters are tested with bumphunter.
#'
#' @param beta beta value matrix.
#' @param pheno phenotype vector with exactly two groups.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param maxClusterGap maximum gap between clusters (bp).
#' @param B number of bumphunter permutations.
#' @param bpSpan bpSpan passed to bumphunter.
#' @param minNum minimum number of probes per block.
#' @param cores number of workers.
#' @param verbose print progress.
#'
#' @return list(Block, clusterInfo, allCLID.v, avbetaCL.m, posCL.m).
#' @export
champ.Block <- function(beta = NULL, pheno = NULL, arraytype = "450K",
                        maxClusterGap = 250000, B = 500L, bpSpan = 250000,
                        minNum = 10L, cores = 3L, verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<< ChAMP.Block START >>>>]")
    .champ_msg("-----------------------------")
    .champ_require("bumphunter", "champ.Block()")
    beta <- .champ_as_beta(beta)
    at <- .champ_arraytype(arraytype)
    pheno <- .champ_pheno(pheno, ncol(beta))
    if (length(unique(pheno)) != 2L) {
        .champ_stop("champ.Block() requires exactly two groups in 'pheno'.")
    }
    if (anyNA(beta)) {
        .champ_msg(sum(is.na(beta)), " NA values detected; imputing.")
        beta <- champ_knn_impute(beta, k = 5L, strict = FALSE, verbose = FALSE)
    }
    pf <- .champ_probe_features(at)
    common <- intersect(rownames(beta), rownames(pf))
    if (length(common) < minNum) .champ_stop("Too few annotated probes.")
    beta <- beta[common, , drop = FALSE]

    RSanno <- pf[common, c("CHR", "MAPINFO", "cgi"), drop = FALSE]
    RSanno <- RSanno[order(suppressWarnings(as.numeric(RSanno$CHR)), RSanno$MAPINFO), ,
                     drop = FALSE]
    sbeta.m <- beta[rownames(RSanno), , drop = FALSE]
    sregion.v <- RSanno$cgi
    .champ_msg("<< Load annotation successfully >>")

    bounded <- function(idx, gap, width) {
        if (!length(idx)) return(integer(0))
        bumphunter::boundedClusterMaker(RSanno$CHR[idx], RSanno$MAPINFO[idx],
                                        assumeSorted = TRUE, maxGap = gap,
                                        maxClusterWidth = width)
    }
    opensea.idx <- which(sregion.v == "opensea")
    island.idx <- which(sregion.v == "island")
    sh.idx <- setdiff(seq_along(sregion.v), c(opensea.idx, island.idx))

    openseaCLID.v <- bounded(opensea.idx, 500, 1500)
    cpgiCLID.v <- bounded(island.idx, 300, 1500)
    shCLID.v <- bounded(sh.idx, 500, 1500)

    allCLID.v <- rep(NA_character_, length(sregion.v))
    allCLID.v[opensea.idx] <- paste0("OS", openseaCLID.v)
    allCLID.v[island.idx] <- paste0("CPGI", cpgiCLID.v)
    allCLID.v[sh.idx] <- paste0("SH", shCLID.v)
    names(allCLID.v) <- rownames(RSanno)

    nf <- length(unique(allCLID.v))
    npg.v <- summary(factor(allCLID.v), maxsum = nf)
    .champ_msg("<< Clusters defined (", nf, ") >>")

    avbetaCL.m <- rowsum(sbeta.m, allCLID.v, reorder = FALSE)
    nf.v <- npg.v[match(rownames(avbetaCL.m), names(npg.v))]
    avbetaCL.m <- avbetaCL.m / as.numeric(nf.v)
    .champ_msg("<< Average beta per cluster calculated >>")

    tmpPOS.v <- RSanno$MAPINFO[match(names(allCLID.v), rownames(RSanno))]
    tmpCHR.v <- .champ_norm_chr(RSanno$CHR[match(names(allCLID.v), rownames(RSanno))])
    tmpCHR.v[tmpCHR.v == "X"] <- "23"
    tmpCHR.v[tmpCHR.v == "Y"] <- "24"
    tmpCHR.v <- suppressWarnings(as.numeric(tmpCHR.v))
    posCL.m <- rowsum(cbind(tmpPOS.v, tmpCHR.v), allCLID.v, reorder = FALSE)
    nf.v <- npg.v[match(rownames(posCL.m), names(npg.v))]
    posCL.m <- posCL.m / as.numeric(nf.v)
    colnames(posCL.m) <- c("Pos", "Chr")

    os.idx <- grep("^OS", rownames(avbetaCL.m))
    if (!length(os.idx)) .champ_stop("No open-sea cluster found; cannot define blocks.")
    avbetaOS.m <- avbetaCL.m[os.idx, , drop = FALSE]
    posOS.m <- posCL.m[os.idx, , drop = FALSE]

    design.m <- data.frame(Intercept = 1,
                           Sample_Group = as.numeric(factor(pheno)) - 1)
    selAUT.idx <- which(posOS.m[, 2L] <= 22)
    if (!length(selAUT.idx)) .champ_stop("No autosomal open-sea cluster found.")
    blocks <- bumphunter::clusterMaker(chr = posOS.m[selAUT.idx, 2L],
                                       pos = posOS.m[selAUT.idx, 1L],
                                       assumeSorted = TRUE, maxGap = maxClusterGap)
    .champ_msg("<< New clusters generated for blocks >>")

    blockPROP.m <- data.frame(
        CHR = stats::aggregate(posOS.m[selAUT.idx, 2L], by = list(blocks),
                               function(x) x[1L])[, 2L],
        stats::aggregate(posOS.m[selAUT.idx, 1L], by = list(blocks),
                         function(x) c(mean(x), min(x), max(x), max(x) - min(x),
                                       length(x)))[, 2L])
    colnames(blockPROP.m) <- c("CHR", "AvPos", "Start", "End", "Size", "NumberOS")
    .champ_msg("<< Block information generated >>")

    bh2.o <- bumphunter::bumphunter(avbetaOS.m[selAUT.idx, , drop = FALSE], design.m,
                                    chr = posOS.m[selAUT.idx, 2L],
                                    pos = posOS.m[selAUT.idx, 1L],
                                    cluster = blocks, coef = 2L, pickCutoff = TRUE,
                                    pickCutoffQ = 0.75, smooth = TRUE,
                                    smoothFunction = bumphunter::loessByCluster,
                                    useWeights = FALSE, B = B, bpSpan = bpSpan,
                                    minNum = minNum)
    .champ_msg("<< Bumphunter finished >>")
    Block <- bh2.o$tab
    if (is.null(Block) || !nrow(Block)) .champ_stop("No differentially methylated block detected.")
    rownames(Block) <- paste0("Block_", seq_len(nrow(Block)))

    .champ_msg("[<<<<< ChAMP.BLOCK END >>>>>]")
    .champ_msg("[===========================]")
    .champ_msg("[You may want to process Block.GUI() next.]\n")
    list(Block = Block, clusterInfo = blockPROP.m, allCLID.v = allCLID.v,
         avbetaCL.m = avbetaCL.m, posCL.m = posCL.m)
}
