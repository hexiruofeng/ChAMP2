test_that("champ.process runs a minimal pipeline and collects results", {
    skip_if_no_champdata()
    beta <- make_beta(300L, 6L)
    pheno <- attr(beta, "group")
    myLoad <- list(beta = beta,
                   pd = data.frame(Sample_Name = colnames(beta),
                                   Sample_Group = pheno,
                                   Slide = rep(c("S1", "S2"), each = 3),
                                   stringsAsFactors = FALSE))
    dir <- tempfile("process")
    res <- suppressMessages(champ.process(
        runload = FALSE, runimpute = FALSE, runQC = FALSE, runnorm = FALSE,
        runSVD = FALSE, runCombat = FALSE, runDMP = TRUE, runDMR = FALSE,
        runBlock = FALSE, runGSEA = FALSE, runEpiMod = FALSE, runCNA = FALSE,
        runRefBase = FALSE, resultsDir = dir, saveStepresults = FALSE,
        arraytype = "450K", verbose = FALSE))
    expect_true("champ.DMP" %in% names(res))
    expect_true(is.list(res[["champ.DMP"]]))
})

test_that("champ.process reports a failing step instead of aborting", {
    skip_if_no_champdata()
    beta <- make_beta(120L, 4L)
    pheno <- c("A", "A", "B", "B")
    myLoad <- list(beta = beta,
                   pd = data.frame(Sample_Name = colnames(beta),
                                   Sample_Group = pheno,
                                   stringsAsFactors = FALSE))
    ## runBlock needs exactly two groups and enough probes; force a failure by
    ## asking for DMRs with minProbes larger than the data
    expect_warning(
        res <- suppressMessages(champ.process(
            runload = FALSE, runimpute = FALSE, runQC = FALSE, runnorm = FALSE,
            runSVD = FALSE, runCombat = FALSE, runDMP = FALSE, runDMR = TRUE,
            runBlock = FALSE, runGSEA = FALSE, runEpiMod = FALSE, runCNA = FALSE,
            runRefBase = FALSE, resultsDir = tempfile(), saveStepresults = FALSE,
            arraytype = "450K", verbose = FALSE)),
        "champ.DMR failed")
    expect_false("champ.DMR" %in% names(res))
})