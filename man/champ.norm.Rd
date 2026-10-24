\name{champ.norm}
\alias{champ.norm}
\title{Normalise methylation data}
\usage{
champ.norm(beta = NULL, rgSet = NULL, mset = NULL, resultsDir =
    "./CHAMP_Normalization/", method = "BMIQ", plotBMIQ = FALSE, arraytype
    = "450K", cores = 3L, verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix.}
  \item{rgSet}{RGChannelSet (required for SWAN and FunctionalNormalization).}
  \item{mset}{MethylSet (required for SWAN).}
  \item{resultsDir}{output directory for BMIQ diagnostic plots.}
  \item{method}{"BMIQ" (default), "PBC", "SWAN" or "FunctionalNormalization".}
  \item{plotBMIQ}{write BMIQ diagnostic PDFs.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{cores}{number of workers for BMIQ.}
  \item{verbose}{print progress.}
}
\description{Normalise methylation data}
\value{normalised beta value matrix.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
