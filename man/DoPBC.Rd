\name{DoPBC}
\alias{DoPBC}
\title{Peak-based correction of type II probe bias}
\usage{
DoPBC(beta.m, design.v, verbose = TRUE)
}
\arguments{
  \item{beta.m}{beta value matrix (probes x samples).}
  \item{design.v}{probe design: 1 = type I, 2 = type II.}
  \item{verbose}{print progress.}
}
\description{Peak-based correction of type II probe bias}
\value{corrected beta value matrix.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
