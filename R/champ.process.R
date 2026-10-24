# =============================================================================
#  champ.process -- the complete ChAMP pipeline in one call
# =============================================================================

#' Run the complete ChAMP pipeline
#'
#' Executes loading, imputation, QC, normalisation, SVD, ComBat, DMP, DMR,
#' block, GSEA, CNA and reference-base correction in the documented order.
#'
#' Unlike the historical implementation, each step is wrapped so that an
#' empty result (for example no significant DMP) does not abort the whole run,
#' and the result list contains exactly the objects that were produced.
#'
#' @param runload,runimpute,runQC,runnorm,runSVD,runCombat,runDMP,runDMR,runBlock,
#'   runGSEA,runEpiMod,runCNA,runRefBase step toggles.
#' @param directory directory with the IDAT files.
#' @param filters which filters to apply (any of "XY", "DetP", "Beads",
#'   "NoCG", "SNP", "MultiHit").
#' @param imputemethod "Combine", "Delete" or "KNN".
#' @param QCplots any of "mdsPlot", "densityPlot", "dendrogram".
#' @param normalizationmethod "BMIQ", "PBC", "SWAN" or "FunctionalNormalization".
#' @param RGEffect include control probe effects in the SVD.
#' @param batchname batch variables corrected by ComBat.
#' @param DMRmethod "Bumphunter", "ProbeLasso" or "DMRcate".
#' @param control,controlGroup CNA controls.
#' @param compare.group two groups to compare.
#' @param adjPVal adjusted p-value cutoff.
#' @param resultsDir output directory.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param PDFplot,Rplot plotting controls.
#' @param cores number of workers.
#' @param saveStepresults save intermediate objects as .rda.
#' @param strict stop on the first error instead of continuing.
#' @param verbose print progress.
#'
#' @return list of the step results.
#' @export
champ.process <- function(runload = TRUE, directory = getwd(),
                          filters = c("XY", "DetP", "Beads", "NoCG", "SNP", "MultiHit"),
                          runimpute = TRUE, imputemethod = "Combine",
                          runQC = TRUE,
                          QCplots = c("mdsPlot", "densityPlot", "dendrogram"),
                          runnorm = TRUE, normalizationmethod = "BMIQ",
                          runSVD = TRUE, RGEffect = FALSE,
                          runCombat = TRUE, batchname = c("Slide"),
                          runDMP = TRUE, runDMR = TRUE, DMRmethod = "Bumphunter",
                          runBlock = TRUE, runGSEA = TRUE, runEpiMod = TRUE,
                          runCNA = TRUE, control = TRUE,
                          controlGroup = "champCtls", runRefBase = FALSE,
                          compare.group = NULL, adjPVal = 0.05,
                          resultsDir = "./CHAMP_RESULT/", arraytype = "450K",
                          PDFplot = TRUE, Rplot = TRUE, cores = 3L,
                          saveStepresults = TRUE, strict = FALSE,
                          verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<< ChAMP.PROCESS START >>>]")
    .champ_msg("-----------------------------")
    at <- .champ_arraytype(arraytype)
    if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)
    .champ_msg("All results will be saved in ", resultsDir)

    RESULT <- list()
    step <- function(name, expr) {
        .champ_msg("\nRunning ", name, "...")
        out <- if (isTRUE(strict)) expr else
            tryCatch(expr, error = function(e) {
                .champ_warn(name, " failed: ", conditionMessage(e))
                structure(list(), class = "champ_step_error",
                          message = conditionMessage(e))
            })
        if (isTRUE(saveStepresults) && !inherits(out, "champ_step_error")) {
            f <- file.path(resultsDir, paste0(name, ".rda"))
            save(out, file = f)
            .champ_msg(name, " result saved in ", f)
        }
        if (!inherits(out, "champ_step_error")) RESULT[[name]] <<- out
        gc()
        .champ_msg("Run ", name, " over.\n")
        invisible(out)
    }

    myLoad <- NULL
    if (isTRUE(runload)) {
        myfilter <- c(XY = "XY", DetP = "DetP", Beads = "Beads", NoCG = "NoCG",
                      SNP = "SNP", MultiHit = "MultiHit") %in% filters
        myLoad <- step("champ.load", champ.load(
            directory = directory, method = "ChAMP", methValue = "B",
            autoimpute = TRUE,
            filterXY = myfilter[["XY"]], filterDetP = myfilter[["DetP"]],
            ProbeCutoff = 0, SampleCutoff = 0.1, detPcut = 0.01,
            filterBeads = myfilter[["Beads"]], beadCutoff = 0.05,
            filterNoCG = myfilter[["NoCG"]], filterSNPs = myfilter[["SNP"]],
            filterMultiHit = myfilter[["MultiHit"]], arraytype = at,
            verbose = verbose))
        if (inherits(myLoad, "champ_step_error")) {
            .champ_stop("champ.load() failed; cannot continue the pipeline.")
        }
    } else {
        if (!exists("myLoad", envir = parent.frame())) {
            .champ_stop("runload = FALSE requires an existing 'myLoad' object.")
        }
        myLoad <- get("myLoad", envir = parent.frame())
    }

    tmppd <- myLoad$pd
    tmpbeta <- myLoad$beta

    myImpute <- NULL
    if (isTRUE(runimpute)) {
        myImpute <- step("champ.impute", champ.impute(
            beta = tmpbeta, pd = tmppd, method = imputemethod, k = 5L,
            ProbeCutoff = 0.2, SampleCutoff = 0.1, verbose = verbose))
        if (!inherits(myImpute, "champ_step_error")) {
            tmpbeta <- myImpute$beta
            tmppd <- myImpute$pd
        }
    }

    if (isTRUE(runQC)) {
        myQC <- c(mdsPlot = "mdsPlot", densityPlot = "densityPlot",
                  dendrogram = "dendrogram") %in% QCplots
        step("champ.QC", champ.QC(
            beta = tmpbeta, pheno = tmppd$Sample_Group,
            mdsPlot = myQC[["mdsPlot"]], densityPlot = myQC[["densityPlot"]],
            dendrogram = myQC[["dendrogram"]], PDFplot = PDFplot, Rplot = Rplot,
            Feature.sel = "None",
            resultsDir = file.path(resultsDir, "CHAMP_QCimages"), verbose = verbose))
    }

    myNorm <- NULL
    if (isTRUE(runnorm)) {
        myNorm <- step("champ.norm", champ.norm(
            beta = tmpbeta, rgSet = myLoad$rgSet, mset = myLoad$mset,
            resultsDir = file.path(resultsDir, "CHAMP_Normalization"),
            method = normalizationmethod, plotBMIQ = PDFplot,
            arraytype = at, cores = cores, verbose = verbose))
        if (!inherits(myNorm, "champ_step_error")) tmpbeta <- myNorm
    }

    if (isTRUE(runSVD)) {
        step("champ.SVD", champ.SVD(
            beta = tmpbeta, rgSet = myLoad$rgSet, pd = tmppd, RGEffect = RGEffect,
            PDFplot = PDFplot, Rplot = Rplot,
            resultsDir = file.path(resultsDir, "CHAMP_SVDimages"), verbose = verbose))
    }

    if (isTRUE(runCombat)) {
        if (!"Slide" %in% colnames(tmppd)) {
            .champ_warn("No 'Slide' column in the sample sheet; ComBat is skipped.")
        } else {
            myCombat <- step("champ.runCombat", champ.runCombat(
                beta = tmpbeta, pd = tmppd, batchname = batchname,
                logitTrans = TRUE, verbose = verbose))
            if (!inherits(myCombat, "champ_step_error")) tmpbeta <- myCombat
        }
    }

    myDMP <- NULL
    if (isTRUE(runDMP)) {
        myDMP <- step("champ.DMP", champ.DMP(
            beta = tmpbeta, pheno = tmppd$Sample_Group, adjPVal = adjPVal,
            adjust.method = "BH", compare.group = compare.group,
            arraytype = at, verbose = verbose))
    }

    myDMR <- NULL
    if (isTRUE(runDMR)) {
        myDMR <- step("champ.DMR", champ.DMR(
            beta = tmpbeta, pheno = tmppd$Sample_Group, compare.group = compare.group,
            arraytype = at, method = DMRmethod, minProbes = 7L,
            adjPvalDmr = adjPVal, cores = cores, maxGap = 300, cutoff = NULL,
            pickCutoff = TRUE, smooth = TRUE, useWeights = FALSE,
            permutations = NULL, B = 250L, nullMethod = "bootstrap",
            meanLassoRadius = 375, minDmrSep = 1000, minDmrSize = 50,
            adjPvalProbe = adjPVal, Rplot = Rplot, PDFplot = PDFplot,
            resultsDir = file.path(resultsDir, "CHAMP_ProbeLasso"),
            rmSNPCH = TRUE, dist = 2, mafcut = 0.05, verbose = verbose))
    }

    myBlock <- NULL
    if (isTRUE(runBlock)) {
        myBlock <- step("champ.Block", champ.Block(
            beta = tmpbeta, pheno = tmppd$Sample_Group, arraytype = at,
            maxClusterGap = 250000, B = 500L, bpSpan = 250000, minNum = 10L,
            cores = cores, verbose = verbose))
    }

    if (isTRUE(runGSEA)) {
        step("champ.GSEA", champ.GSEA(
            beta = tmpbeta,
            DMP = if (inherits(myDMP, "champ_step_error")) NULL else myDMP,
            DMR = if (inherits(myDMR, "champ_step_error")) NULL else myDMR,
            CpGlist = NULL, Genelist = NULL, arraytype = at, adjPval = adjPVal,
            verbose = verbose))
    }

    if (isTRUE(runEpiMod)) {
        step("champ.ebGSEA", champ.ebGSEA(
            beta = tmpbeta, pheno = tmppd$Sample_Group, adjPval = adjPVal,
            arraytype = at, cores = 1L, verbose = verbose))
    }

    if (isTRUE(runCNA)) {
        step("champ.CNA", champ.CNA(
            intensity = myLoad$intensity, pheno = myLoad$pd$Sample_Group,
            resultsDir = file.path(resultsDir, "CHAMP_CNA"), control = control,
            controlGroup = controlGroup, sampleCNA = TRUE, groupFreqPlots = TRUE,
            Rplot = Rplot, PDFplot = PDFplot, freqThreshold = 0.3,
            arraytype = at, verbose = verbose))
    }

    if (isTRUE(runRefBase)) {
        step("champ.refbase", champ.refbase(beta = tmpbeta, arraytype = at,
                                            verbose = verbose))
    }

    .champ_msg("[<<<< ChAMP.PROCESS END >>>>]")
    .champ_msg("[===========================]")
    RESULT
}
