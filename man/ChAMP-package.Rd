\name{ChAMP-package}
\alias{ChAMP-package}
\alias{ChAMP}
\docType{package}
\title{ChAMP2: Chip Analysis Methylation Pipeline}
\description{
  Quality control, normalisation, differential methylation, copy number and
  gene set enrichment analysis for Illumina HumanMethylation 450K, EPIC and
  EPICv2 arrays.  This fork (ChAMP2) installs and runs on R 4.2.3 to 4.6.2.
}
\details{
  The main entry point is \code{\link{champ.process}}; individual steps can be
  run with \code{\link{champ.load}}, \code{\link{champ.filter}},
  \code{\link{champ.QC}}, \code{\link{champ.norm}}, \code{\link{champ.SVD}},
  \code{\link{champ.runCombat}}, \code{\link{champ.DMP}}, \code{\link{champ.DMR}},
  \code{\link{champ.Block}}, \code{\link{champ.GSEA}}, \code{\link{champ.ebGSEA}},
  \code{\link{champ.CNA}} and \code{\link{champ.refbase}}.
  Run \code{\link{champ_check_env}} to see which optional features are available
  and \code{\link{champ_install_deps}} to install the missing packages.
}
\keyword{package}
