\name{champ.ComBat}
\alias{champ.ComBat}
\title{ComBat batch correction}
\usage{
champ.ComBat(dat, batch, mod, numCovs = NULL, par.prior = TRUE, prior.plots
    = FALSE, mean.only = FALSE)
}
\arguments{
  \item{dat}{data matrix (features x samples).}
  \item{batch}{batch vector.}
  \item{mod}{model matrix of covariates (may include a "Batch" column, which is removed because ComBat adds the batch term itself).}
  \item{numCovs}{unused, kept for backwards compatibility.}
  \item{par.prior}{use the parametric empirical Bayes prior.}
  \item{prior.plots}{diagnostic plots.}
  \item{mean.only}{only correct the mean (sva >= 3.30).}
}
\description{ComBat batch correction}
\value{the corrected matrix.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
