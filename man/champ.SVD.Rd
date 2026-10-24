\name{champ.SVD}
\alias{champ.SVD}
\title{SVD analysis of phenotype associations}
\usage{
champ.SVD(beta = NULL, rgSet = NULL, pd = NULL, RGEffect = FALSE, PDFplot =
    TRUE, Rplot = TRUE, resultsDir = "./CHAMP_SVDimages/", verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix.}
  \item{rgSet}{optional RGChannelSet used to add control probe effects.}
  \item{pd}{sample sheet.}
  \item{RGEffect}{include control probe intensities as covariates.}
  \item{PDFplot}{Write plots to PDF files in resultsDir.}
  \item{Rplot}{Draw plots on the current device.}
  \item{resultsDir}{output directory.}
  \item{verbose}{print progress.}
}
\description{Computes the singular value decomposition of the centred beta matrix, estimates the number of informative components with random matrix theory and tests each component against every factor in \\code\{pd\} (Kruskal-Wallis for categorical variables, linear regression for numeric ones).}
\value{matrix of p-values (components x factors).}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
