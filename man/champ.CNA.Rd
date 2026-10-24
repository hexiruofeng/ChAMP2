\name{champ.CNA}
\alias{champ.CNA}
\title{Copy number analysis}
\usage{
champ.CNA(intensity = NULL, pheno = NULL, control = TRUE, controlGroup =
    "champCtls", sampleCNA = TRUE, groupFreqPlots = TRUE, Rplot = FALSE,
    PDFplot = TRUE, freqThreshold = 0.3, resultsDir = "./CHAMP_CNA",
    arraytype = "450K", verbose = TRUE)
}
\arguments{
  \item{intensity}{total intensity matrix (probes x samples).}
  \item{pheno}{phenotype vector.}
  \item{control}{use control samples for the log ratio.}
  \item{controlGroup}{name of the control group ("champCtls" uses the bundled blood controls).}
  \item{sampleCNA}{segment every sample.}
  \item{groupFreqPlots}{draw per-group gain/loss frequency plots.}
  \item{Rplot}{Draw plots on the current device.}
  \item{PDFplot}{Write plots to PDF files in resultsDir.}
  \item{freqThreshold}{frequency threshold for the group plots.}
  \item{resultsDir}{output directory.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{verbose}{print progress.}
}
\description{Quantile-normalises total intensities, computes log2 ratios against control samples (ChAMP blood controls by default) and segments each sample with DNAcopy.}
\value{list(sampleResult, groupResult).}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
