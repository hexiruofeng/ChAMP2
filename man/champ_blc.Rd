\name{champ_blc}
\alias{champ_blc}
\title{Beta latent class model (EM)}
\usage{
champ_blc(Y, w, maxiter = 25L, tol = 1e-06, weights = NULL, verbose = FALSE)
}
\arguments{
  \item{Y}{numeric matrix (n x J), typically a single column of beta values.}
  \item{w}{numeric matrix (n x K) of initial posterior class memberships.}
  \item{maxiter}{maximum EM iterations.}
  \item{tol}{convergence tolerance on mu.}
  \item{weights}{case weights.}
  \item{verbose}{print the convergence criterion at each iteration.}
}
\description{Re-implementation of the beta mixture EM used by BMIQ.  Returns the same components as \\code\{RPMM::blc\}: \\code\{a\}, \\code\{b\}, \\code\{eta\}, \\code\{mu\}, \\code\{w\} and the log-likelihood.}
\value{list with a, b (K x J), eta (K), mu (K x J), w (n x K), llike.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
