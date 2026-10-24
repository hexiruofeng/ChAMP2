\name{QC.GUI}
\alias{QC.GUI}
\title{Interactive quality control}
\usage{
QC.GUI(beta = NULL, pheno = NULL, arraytype = "450K", verbose = TRUE,
    runApp = TRUE, numPositions = 1000L)
}
\arguments{
  \item{beta}{beta value matrix (probes x samples).}
  \item{pheno}{phenotype vector, one entry per column of \\code\{beta\}.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{verbose}{print progress messages.}
  \item{runApp}{launch the application; FALSE returns the shiny.appobj.}
  \item{numPositions}{number of most variable probes used for the MDS/PCA.}
}
\description{Interactive quality control}
\value{invisibly, the shiny.appobj returned by \\code\{shiny::runApp\}.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
