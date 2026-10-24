\name{champ.read}
\alias{champ.read}
\title{Read a beta matrix and a sample sheet from delimited files}
\usage{
champ.read(betaFile = "beta.txt", sampleSheet = "sampleSheet.txt",
    resultsDir = NULL, sep = "\\t", verbose = TRUE)
}
\arguments{
  \item{betaFile}{tab separated beta matrix file.}
  \item{sampleSheet}{tab separated sample sheet.}
  \item{resultsDir}{output directory (created when missing).}
  \item{sep}{field separator.}
  \item{verbose}{print progress.}
}
\description{Read a beta matrix and a sample sheet from delimited files}
\value{list(beta, pd).}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
