\name{champ.ebGSEA}
\alias{champ.ebGSEA}
\title{Empirical Bayes gene set enrichment}
\usage{
champ.ebGSEA(beta = NULL, pheno = NULL, minN = 5L, adjPval = 0.05,
    arraytype = "450K", cores = 1L, verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix.}
  \item{pheno}{phenotype vector (binary factor or numeric covariate).}
  \item{minN}{minimum number of genes per set.}
  \item{adjPval}{adjusted p-value cutoff.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{cores}{number of workers.}
  \item{verbose}{print progress.}
}
\description{Empirical Bayes gene set enrichment}
\value{list with GSEA, EnrichGene and gtResult.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
