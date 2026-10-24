\name{champ_install_deps}
\alias{champ_install_deps}
\title{Install the packages ChAMP needs}
\usage{
champ_install_deps(optional = TRUE, binary = TRUE, upgrade = "never", force
    = FALSE, verbose = TRUE)
}
\arguments{
  \item{optional}{also install the Suggests packages.}
  \item{binary}{use binary packages when available (much faster on Windows and macOS).}
  \item{upgrade}{"never" (default), "default" or "always".}
  \item{force}{reinstall packages that are already present.}
  \item{verbose}{print progress.}
}
\description{Installs every missing required package (Bioconductor first, then CRAN) and, unless \\code\{optional = FALSE\}, the suggested packages that enable the optional features (minfi loading, DMRcate backend, RPMM/impute/kpmt legacy backends, the Shiny GUIs, ...). The function is deliberately forgiving: a package that cannot be installed (for example because the current Bioconductor release dropped it) produces a warning and the installation continues, because ChAMP bundles a fallback for each of those algorithms.}
\value{invisibly, a list with installed, skipped and failed packages.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
