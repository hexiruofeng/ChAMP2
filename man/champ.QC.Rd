\name{champ.QC}
\alias{champ.QC}
\title{Quality control plots}
\usage{
champ.QC(beta = NULL, pheno = NULL, mdsPlot = TRUE, densityPlot = TRUE,
    dendrogram = TRUE, PDFplot = TRUE, Rplot = TRUE, Feature.sel = "None",
    resultsDir = "./CHAMP_QCimages/", verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix.}
  \item{pheno}{phenotype vector.}
  \item{mdsPlot}{See the package vignette for details.}
  \item{densityPlot}{See the package vignette for details.}
  \item{dendrogram}{See the package vignette for details.}
  \item{PDFplot}{write PDF files into resultsDir.}
  \item{Rplot}{draw on the current device.}
  \item{Feature.sel}{"None" or "SVD" probe selection for the dendrogram.}
  \item{resultsDir}{output directory.}
  \item{verbose}{print progress.}
}
\description{Produces MDS, density and dendrogram plots before and/or after normalisation.  The plotting helpers are self-contained (no minfi needed).}
\value{invisibly, the hclust object when a dendrogram is produced.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
