\name{champ.DMP}
\alias{champ.DMP}
\title{Differential methylation analysis at the probe level}
\usage{
champ.DMP(beta = NULL, pheno = NULL, compare.group = NULL, adjPVal = 0.05,
    adjust.method = "BH", arraytype = "450K", verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix (probes x samples).}
  \item{pheno}{phenotype vector (character, factor or numeric covariate).}
  \item{compare.group}{two group labels to compare (optional).}
  \item{adjPVal}{adjusted p-value cutoff used for reporting.}
  \item{adjust.method}{multiple testing correction ("BH" by default).}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{verbose}{print progress.}
}
\description{For a categorical phenotype every pair of groups is compared with limma (moderated t-test); for a numeric phenotype a linear regression is fitted. Results are annotated with the ChAMP probe annotation.}
\value{named list of data.frames, one per comparison (or "NumericVariable").}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
