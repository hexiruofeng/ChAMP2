\name{champ.impute}
\alias{champ.impute}
\title{Impute missing beta values}
\usage{
champ.impute(beta = NULL, pd = NULL, method = "Combine", k = 5L,
    ProbeCutoff = 0.2, SampleCutoff = 0.1, verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix (probes x samples).}
  \item{pd}{sample sheet; rows are subset alongside the matrix.}
  \item{method}{"Combine", "Delete" or "KNN".}
  \item{k}{number of neighbours for KNN.}
  \item{ProbeCutoff}{maximum fraction of missing values per probe (Combine).}
  \item{SampleCutoff}{maximum fraction of missing values per sample (Combine).}
  \item{verbose}{print progress.}
}
\description{Three strategies are available: \\describe\{ \\item\{Combine\}\{remove samples and probes with too many missing values, then KNN-impute the rest (recommended for large data sets).\} \\item\{Delete\}\{remove any probe with at least one missing value (recommended for small data sets).\} \\item\{KNN\}\{KNN-impute every missing value.\} \}}
\value{list(beta, pd) when pd is given, otherwise the imputed matrix.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
