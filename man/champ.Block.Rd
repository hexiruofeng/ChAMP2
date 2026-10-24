\name{champ.Block}
\alias{champ.Block}
\title{Differentially methylated blocks}
\usage{
champ.Block(beta = NULL, pheno = NULL, arraytype = "450K", maxClusterGap =
    250000, B = 500L, bpSpan = 250000, minNum = 10L, cores = 3L, verbose =
    TRUE)
}
\arguments{
  \item{beta}{beta value matrix.}
  \item{pheno}{phenotype vector with exactly two groups.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{maxClusterGap}{maximum gap between clusters (bp).}
  \item{B}{number of bumphunter permutations.}
  \item{bpSpan}{bpSpan passed to bumphunter.}
  \item{minNum}{minimum number of probes per block.}
  \item{cores}{number of workers.}
  \item{verbose}{print progress.}
}
\description{Probes are aggregated into CpG-island / open-sea / shore-shelf clusters and blocks of spatially adjacent clusters are tested with bumphunter.}
\value{list(Block, clusterInfo, allCLID.v, avbetaCL.m, posCL.m).}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
