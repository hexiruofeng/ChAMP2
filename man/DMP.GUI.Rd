\name{DMP.GUI}
\alias{DMP.GUI}
\title{Explore differentially methylated positions in a Shiny application}
\usage{
DMP.GUI(DMP = NULL, beta = NULL, pheno = NULL, cutgroupnumber = 4,
    arraytype = "450K", verbose = TRUE, runApp = TRUE)
}
\arguments{
  \item{DMP}{data.frame (or list of data.frames) returned by champ.DMP().}
  \item{beta}{beta value matrix used for the DMP test (probes x samples).}
  \item{pheno}{phenotype vector, one entry per column of \\code\{beta\}. A numeric phenotype is binned into \\code\{cutgroupnumber\} groups.}
  \item{cutgroupnumber}{number of groups used when \\code\{pheno\} is numeric.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2"; used for missing annotation.}
  \item{verbose}{print progress messages.}
  \item{runApp}{launch the application; when FALSE the shiny.appobj is returned instead (useful for tests).}
}
\description{Opens a Shiny/plotly interface with the DMP table, heatmap, feature/CpG composition, per-gene profile and per-CpG distribution of a comparison produced by \\code\{\\link\{champ.DMP\}\}.}
\value{invisibly, the shiny.appobj returned by \\code\{shiny::runApp\}.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
