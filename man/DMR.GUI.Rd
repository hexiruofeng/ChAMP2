\name{DMR.GUI}
\alias{DMR.GUI}
\title{Explore differentially methylated regions in a Shiny application}
\usage{
DMR.GUI(DMR = NULL, beta = NULL, pheno = NULL, arraytype = "450K", verbose
    = TRUE, runApp = TRUE)
}
\arguments{
  \item{DMR}{data.frame (or list) returned by champ.DMR().}
  \item{beta}{beta value matrix used for the DMR analysis.}
  \item{pheno}{phenotype vector, one entry per column of \\code\{beta\}.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{verbose}{print progress messages.}
  \item{runApp}{launch the application; FALSE returns the shiny.appobj.}
}
\description{Explore differentially methylated regions in a Shiny application}
\value{invisibly, the shiny.appobj returned by \\code\{shiny::runApp\}.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
