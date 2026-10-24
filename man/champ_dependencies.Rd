\name{champ_dependencies}
\alias{champ_dependencies}
\title{ChAMP dependency manifest}
\usage{
champ_dependencies(which = c("all", "Imports", "Suggests"))
}
\arguments{
  \item{which}{"Imports", "Suggests" or "all".}
}
\description{ChAMP dependency manifest}
\value{data.frame with package, role, source, min_version and feature.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
