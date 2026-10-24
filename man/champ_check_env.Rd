\name{champ_check_env}
\alias{champ_check_env}
\title{Report the status of every ChAMP dependency}
\usage{
champ_check_env(quiet = FALSE)
}
\arguments{
  \item{quiet}{do not print the report.}
}
\description{Prints (and invisibly returns) a table with the installed version of every package ChAMP can use, whether it is required, and which feature it enables.}
\value{invisibly, a data.frame.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
