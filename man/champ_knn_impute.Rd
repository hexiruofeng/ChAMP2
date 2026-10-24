\name{champ_knn_impute}
\alias{champ_knn_impute}
\title{K-Nearest-Neighbour imputation}
\usage{
champ_knn_impute(data, k = 5L, rowmax = 0.5, colmax = 0.8, maxp = 1500L,
    use_impute = TRUE, strict = TRUE, verbose = FALSE)
}
\arguments{
  \item{data}{numeric matrix, rows = features, columns = samples.}
  \item{k}{number of neighbours.}
  \item{rowmax}{features with more than this fraction of missing values are filled with the column mean (as in \\code\{impute.knn\}).}
  \item{colmax}{samples with more than this fraction of missing values are rejected; for compatibility with the historical behaviour this is an error, unless \\code\{strict = FALSE\}.}
  \item{maxp}{maximum number of rows used to estimate sample distances (\\code\{impute.knn\} splits the matrix instead; results are equivalent for the purpose of ChAMP).}
  \item{use_impute}{use the \\code\{impute\} package when installed (default).}
  \item{strict}{error out on samples exceeding \\code\{colmax\} (default TRUE, matching \\code\{impute.knn\}); when FALSE such samples are mean-filled.}
  \item{verbose}{print progress messages.}
}
\description{Drop-in replacement for \\code\{impute::impute.knn\}.  Rows are features (probes) and columns are samples; missing entries of a feature are imputed from the feature's values in the \\code\{k\} most similar samples. The original \\code\{impute\} package is used when available, so that existing pipelines keep bit-identical results.  Otherwise a vectorised pure-R implementation with the same arguments is used.}
\value{the imputed matrix, with dimnames preserved.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
