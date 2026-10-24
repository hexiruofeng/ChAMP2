\name{champ_kpmt}
\alias{champ_kpmt}
\title{Known Population Median Test}
\usage{
champ_kpmt(pop, obs = NULL, med = NULL, size = NULL, tail = "lower")
}
\arguments{
  \item{pop}{named numeric vector (or one-column matrix) of population values.}
  \item{obs}{character vector of population member names that were observed.}
  \item{med}{pre-computed minimal median.}
  \item{size}{sample size (only needed when \\code\{med\} is supplied).}
  \item{tail}{one of "lower", "upper", "two-sided".}
}
\description{Exact finite-population test for whether the minimal median of an observed subset is shifted relative to the whole population.  This is a re-implementation of \\code\{kpmt::kpmt\} (CRAN, archived 2025) and returns a data frame with the same columns, so existing code that indexes the fourth column (\\code\{p.value\}) keeps working. The p-value is the exact combinatorial probability \\deqn\{P(\\tilde\{X\} \\le w) = \\frac\{s!(N-s)!\}\{N!(m-1)!(s-m)!\} \\sum_\{x \\in U\} \\frac\{(x-1)!(N-x)!\}\{(x-m)!(N-x-s+m)!\}\} with \\eqn\{m = \\lceil s/2 \\rceil\} and \\eqn\{U\} the population ranks that are below the observed minimal median and still leave room for the sample.}
\value{data.frame with columns feature, sample, logp, p.value, median.sample, median.all, median.other.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
