\name{champ.import}
\alias{champ.import}
\title{Import Illumina IDAT files using the ChAMP annotation}
\usage{
champ.import(directory = getwd(), offset = 100, arraytype = "450K", verbose
    = TRUE)
}
\arguments{
  \item{directory}{directory containing the sample sheet (csv) and IDAT files.}
  \item{offset}{offset added to the denominator when computing beta values.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2" ("EPIC" is an alias of "EPICv2").}
  \item{verbose}{print progress.}
}
\description{Reads a sample sheet and the matching Grn/Red IDAT files, and returns beta values, M values, intensities, detection p-values and bead counts.}
\value{list with beta, M, pd, intensity, detP, beadcount, Meth, UnMeth, arraytype.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
