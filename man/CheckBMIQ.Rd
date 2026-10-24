\name{CheckBMIQ}
\alias{CheckBMIQ}
\title{Check a BMIQ normalised profile}
\usage{
CheckBMIQ(beta.v, design.v, pnbeta.v)
}
\arguments{
  \item{beta.v}{original beta values.}
  \item{design.v}{probe design (1/2).}
  \item{pnbeta.v}{BMIQ normalised beta values.}
}
\description{Check a BMIQ normalised profile}
\value{invisibly, a data.frame with the observed and normalised densities.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
