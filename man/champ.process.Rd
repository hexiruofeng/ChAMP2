\name{champ.process}
\alias{champ.process}
\title{Run the complete ChAMP pipeline}
\usage{
champ.process(runload = TRUE, directory = getwd(), filters = c("XY",
    "DetP", "Beads", "NoCG", "SNP", "MultiHit"), runimpute = TRUE,
    imputemethod = "Combine", runQC = TRUE, QCplots = c("mdsPlot",
    "densityPlot", "dendrogram"), runnorm = TRUE, normalizationmethod =
    "BMIQ", runSVD = TRUE, RGEffect = FALSE, runCombat = TRUE, batchname =
    c("Slide"), runDMP = TRUE, runDMR = TRUE, DMRmethod = "Bumphunter",
    runBlock = TRUE, runGSEA = TRUE, runEpiMod = TRUE, runCNA = TRUE,
    control = TRUE, controlGroup = "champCtls", runRefBase = FALSE,
    compare.group = NULL, adjPVal = 0.05, resultsDir = "./CHAMP_RESULT/",
    arraytype = "450K", PDFplot = TRUE, Rplot = TRUE, cores = 3L,
    saveStepresults = TRUE, strict = FALSE, verbose = TRUE)
}
\arguments{
  \item{runload}{See the package vignette for details.}
  \item{directory}{directory with the IDAT files.}
  \item{filters}{which filters to apply (any of "XY", "DetP", "Beads", "NoCG", "SNP", "MultiHit").}
  \item{runimpute}{See the package vignette for details.}
  \item{imputemethod}{"Combine", "Delete" or "KNN".}
  \item{runQC}{See the package vignette for details.}
  \item{QCplots}{any of "mdsPlot", "densityPlot", "dendrogram".}
  \item{runnorm}{See the package vignette for details.}
  \item{normalizationmethod}{"BMIQ", "PBC", "SWAN" or "FunctionalNormalization".}
  \item{runSVD}{See the package vignette for details.}
  \item{RGEffect}{include control probe effects in the SVD.}
  \item{runCombat}{See the package vignette for details.}
  \item{batchname}{batch variables corrected by ComBat.}
  \item{runDMP}{See the package vignette for details.}
  \item{runDMR}{See the package vignette for details.}
  \item{DMRmethod}{"Bumphunter", "ProbeLasso" or "DMRcate".}
  \item{runBlock}{See the package vignette for details.}
  \item{runGSEA}{See the package vignette for details.}
  \item{runEpiMod}{See the package vignette for details.}
  \item{runCNA}{See the package vignette for details.}
  \item{control}{Use control samples to compute the CNA log ratios.}
  \item{controlGroup}{Name of the control group; 'champCtls' uses the bundled blood controls.}
  \item{runRefBase}{See the package vignette for details.}
  \item{compare.group}{two groups to compare.}
  \item{adjPVal}{adjusted p-value cutoff.}
  \item{resultsDir}{output directory.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{PDFplot}{Write plots to PDF files in resultsDir.}
  \item{Rplot}{Draw plots on the current device.}
  \item{cores}{number of workers.}
  \item{saveStepresults}{save intermediate objects as .rda.}
  \item{strict}{stop on the first error instead of continuing.}
  \item{verbose}{print progress.}
}
\description{Executes loading, imputation, QC, normalisation, SVD, ComBat, DMP, DMR, block, GSEA, CNA and reference-base correction in the documented order. Unlike the historical implementation, each step is wrapped so that an empty result (for example no significant DMP) does not abort the whole run, and the result list contains exactly the objects that were produced.}
\value{list of the step results.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
