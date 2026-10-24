\name{champ.refbase}
\alias{champ.refbase}
\title{Estimate and correct blood cell proportions}
\usage{
champ.refbase(beta = NULL, arraytype = "450K", verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix.}
  \item{arraytype}{"450K", "EPICv1" (both use the 450K reference) or "27K".}
  \item{verbose}{print progress.}
}
\description{Uses a reference-based (constrained least squares) deconvolution with the Reinius et al. (450K) or Accomando et al. (27K) reference profiles and regresses out the estimated proportions.}
\value{list(CorrectedBeta, CellFraction).}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
