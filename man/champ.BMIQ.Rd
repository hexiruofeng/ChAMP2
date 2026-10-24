\name{champ.BMIQ}
\alias{champ.BMIQ}
\title{BMIQ normalisation of one sample}
\usage{
champ.BMIQ(beta.v, design.v, nL = 3, doH = TRUE, nfit = 10000, th1.v =
    c(0.2, 0.75), th2.v = NULL, niter = 5, tol = 0.001, plots = FALSE,
    plotDir = ".", sampleID = 1, verbose = TRUE)
}
\arguments{
  \item{beta.v}{beta values of one sample (no NA).}
  \item{design.v}{probe design: 1 = type I, 2 = type II.}
  \item{nL}{number of mixture components (BMIQ requires 3).}
  \item{doH}{correct hemimethylated type II probes.}
  \item{nfit}{number of probes used to fit the mixture (capped at the number of available probes).}
  \item{th1.v}{initial thresholds for type I probes.}
  \item{th2.v}{initial thresholds for type II probes (estimated when NULL).}
  \item{niter}{number of EM iterations.}
  \item{tol}{EM tolerance.}
  \item{plots}{write diagnostic PDFs.}
  \item{plotDir}{directory for the diagnostic PDFs.}
  \item{sampleID}{sample label used in file names.}
  \item{verbose}{print progress.}
}
\description{BMIQ normalisation of one sample}
\value{list with nbeta (normalised beta), class1, class2, av1, av2, hf, th1, th2.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
