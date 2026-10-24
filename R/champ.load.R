# =============================================================================
#  champ.load -- load an Illumina methylation experiment
#
#  Two back-ends:
#    * "ChAMP"  : champ.import() + champ.filter()  (no minfi needed)
#    * "minfi"  : minfi::read.metharray.exp() + preprocessRaw(), then the same
#                 filtering engine, so both paths behave identically.
# =============================================================================

#' Bead counts from an RGChannelSet
#' @keywords internal
.champ_beadcount <- function(rgSet) {
    nb <- minfi::getNBeads(rgSet)
    loci <- minfi::getManifestInfo(rgSet, "locusNames")
    bc <- matrix(NA_real_, nrow = length(loci), ncol = ncol(rgSet),
                 dimnames = list(loci, minfi::sampleNames(rgSet)))
    typeII <- minfi::getProbeInfo(rgSet, type = "II")
    typeI  <- minfi::getProbeInfo(rgSet, type = "I")
    addrA <- typeII$AddressA
    okII <- addrA %in% rownames(nb)
    bc[typeII$Name[okII], ] <- nb[addrA[okII], , drop = FALSE]

    bcA <- bcB <- bc
    okA <- typeI$AddressA %in% rownames(nb)
    okB <- typeI$AddressB %in% rownames(nb)
    bcA[typeI$Name[okA], ] <- nb[typeI$AddressA[okA], , drop = FALSE]
    bcB[typeI$Name[okB], ] <- nb[typeI$AddressB[okB], , drop = FALSE]
    bad <- (bcA < 3) | (bcB < 3)
    bad[is.na(bad)] <- TRUE
    bc[bad] <- NA_real_
    bc
}

#' Load methylation array data
#'
#' @param directory directory containing the IDAT files (and, for method
#'   "minfi", a sample sheet).
#' @param method "ChAMP" (annotation based, default) or "minfi".
#' @param methValue "B" (beta) or "M" (M-value).
#' @param autoimpute impute remaining missing values.
#' @param filterDetP,filterBeads,filterNoCG,filterSNPs,filterMultiHit,filterXY
#'   filtering toggles, see champ.filter().
#' @param ProbeCutoff,SampleCutoff,detPcut,beadCutoff filtering thresholds.
#' @param population SNP mask population.
#' @param force passed to minfi::read.metharray.exp().
#' @param arraytype "450K", "EPICv1" or "EPICv2" ("EPIC" = "EPICv2").
#' @param verbose print progress.
#'
#' @return list with beta (or M), pd, intensity, detP and, for method
#'   "minfi", mset and rgSet.
#' @export
champ.load <- function(directory = getwd(), method = "ChAMP", methValue = "B",
                       autoimpute = TRUE, filterDetP = TRUE, ProbeCutoff = 0,
                       SampleCutoff = 0.1, detPcut = 0.01, filterBeads = TRUE,
                       beadCutoff = 0.05, filterNoCG = TRUE, filterSNPs = TRUE,
                       population = NULL, filterMultiHit = TRUE, filterXY = TRUE,
                       force = FALSE, arraytype = "450K", verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<< ChAMP.LOAD START >>>>>]")
    .champ_msg("-----------------------------")
    at <- .champ_arraytype(arraytype)
    method <- match.arg(method, c("ChAMP", "minfi"))

    if (identical(method, "minfi")) {
        .champ_require("minfi", "champ.load(method = 'minfi')")
        .champ_msg("\n[ Loading Data with Minfi Method ]")
        .champ_msg("----------------------------------")
        .champ_msg("Loading data from ", directory)

        targets <- .champ_try(suppressWarnings(minfi::read.metharray.sheet(directory)))
        if (is.null(targets) || !nrow(targets)) {
            .champ_stop("minfi could not read a sample sheet in ", directory,
                        ".\n  Provide a sample sheet (csv) or use method = 'ChAMP'.")
        }
        rgSet <- minfi::read.metharray.exp(targets = targets, extended = TRUE,
                                           force = force)
        if (identical(at, "EPICv1")) {
            rgSet@annotation <- c(array = "IlluminaHumanMethylationEPIC",
                                  annotation = "ilm10b4.hg19")
        } else if (identical(at, "EPICv2")) {
            rgSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2",
                                  annotation = "20a1.hg38")
        }
        if ("Sample_Name" %in% colnames(targets)) {
            minfi::sampleNames(rgSet) <- as.character(targets$Sample_Name)
        }
        pd <- as.data.frame(minfi::pData(rgSet), stringsAsFactors = FALSE)
        if (!"Sample_Name" %in% colnames(pd)) pd$Sample_Name <- colnames(rgSet)
        mset <- minfi::preprocessRaw(rgSet)
        detP <- minfi::detectionP(rgSet)
        .champ_msg("<< Read DataSet Success. >>\n")

        beta <- if (identical(methValue, "B")) {
            minfi::getBeta(mset, "Illumina")
        } else {
            minfi::getM(mset)
        }
        beta[detP >= detPcut] <- NA_real_

        numfail <- matrix(colMeans(is.na(beta)))
        rownames(numfail) <- colnames(detP)
        colnames(numfail) <- "Failed CpG Fraction."
        print(numfail)

        beadcount <- if (isTRUE(filterBeads)) {
            .champ_try(.champ_beadcount(rgSet))
        } else NULL

        res <- champ.filter(beta = beta, pd = pd, detP = detP,
                            beadcount = beadcount, autoimpute = autoimpute,
                            filterDetP = filterDetP, ProbeCutoff = ProbeCutoff,
                            SampleCutoff = SampleCutoff, detPcut = detPcut,
                            filterBeads = filterBeads && !is.null(beadcount),
                            beadCutoff = beadCutoff, filterNoCG = filterNoCG,
                            filterSNPs = filterSNPs, population = population,
                            filterMultiHit = filterMultiHit, filterXY = filterXY,
                            arraytype = at, verbose = verbose)

        samples <- colnames(res$beta)
        probes  <- rownames(res$beta)
        rgSet <- rgSet[, samples]
        mset <- mset[minfi::featureNames(mset) %in% probes, ]
        mset <- mset[match(probes, minfi::featureNames(mset)), ]
        intensity <- minfi::getMeth(mset) + minfi::getUnmeth(mset)
        detP <- detP[probes, samples, drop = FALSE]

        .champ_msg("The analysis will proceed with ", nrow(res$beta), " probes and ",
                   ncol(res$beta), " samples.\n")
        .champ_msg("[<<<<< ChAMP.LOAD END >>>>>>]")
        .champ_msg("[===========================]")
        .champ_msg("[You may want to process champ.QC() next.]\n")
        return(list(mset = mset, rgSet = rgSet, pd = res$pd,
                    intensity = intensity, beta = res$beta, detP = detP,
                    arraytype = at))
    }

    .champ_msg("\n[ Loading Data with ChAMP Method ]")
    .champ_msg("----------------------------------")
    .champ_msg("Note that the ChAMP method does NOT return rgSet or mset, so SWAN ",
               "and FunctionalNormalization are unavailable; BMIQ and PBC work fine.\n")

    myImport <- champ.import(directory, arraytype = at, verbose = verbose)
    res <- champ.filter(
        beta = if (identical(methValue, "B")) myImport$beta else NULL,
        M    = if (identical(methValue, "B")) NULL else myImport$M,
        pd = myImport$pd, intensity = myImport$intensity,
        Meth = myImport$Meth, UnMeth = myImport$UnMeth,
        detP = myImport$detP, beadcount = myImport$beadcount,
        autoimpute = autoimpute, filterDetP = filterDetP,
        ProbeCutoff = ProbeCutoff, SampleCutoff = SampleCutoff,
        detPcut = detPcut, filterBeads = filterBeads,
        beadCutoff = beadCutoff, filterNoCG = filterNoCG,
        filterSNPs = filterSNPs, population = population,
        filterMultiHit = filterMultiHit, filterXY = filterXY,
        arraytype = at, verbose = verbose)

    .champ_msg("[<<<<< ChAMP.LOAD END >>>>>>]")
    .champ_msg("[===========================]")
    .champ_msg("[You may want to process champ.QC() next.]\n")
    res$arraytype <- at
    res
}
