\name{Block.GUI}
\alias{Block.GUI}
\title{Explore differentially methylated blocks in a Shiny application}
\usage{
Block.GUI(Block = NULL, beta = NULL, pheno = NULL, arraytype = "450K",
    verbose = TRUE, runApp = TRUE)
}
\arguments{
  \item{Block}{data.frame (or the list returned by champ.Block()).}
  \item{beta}{beta value matrix used for the block analysis.}
  \item{pheno}{phenotype vector, one entry per column of \\code\{beta\}.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{verbose}{print progress messages.}
  \item{runApp}{launch the application; FALSE returns the shiny.appobj.}
}
\description{Explore differentially methylated blocks in a Shiny application}
\value{invisibly, the shiny.appobj returned by \\code\{shiny::runApp\}.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
