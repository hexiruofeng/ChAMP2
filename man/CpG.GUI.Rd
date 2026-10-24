\name{CpG.GUI}
\alias{CpG.GUI}
\title{Explore the annotation of a CpG list}
\usage{
CpG.GUI(CpG = NULL, arraytype = "450K", verbose = TRUE, runApp = TRUE)
}
\arguments{
  \item{CpG}{character vector of CpG identifiers. NULL uses every probe of the array.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{verbose}{print progress messages.}
  \item{runApp}{launch the application; FALSE returns the shiny.appobj.}
}
\description{Explore the annotation of a CpG list}
\value{invisibly, the shiny.appobj returned by \\code\{shiny::runApp\}.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
