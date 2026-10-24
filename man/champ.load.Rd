\name{champ.load}
\alias{champ.load}
\title{Load methylation array data}
\usage{
champ.load(directory = getwd(), method = "ChAMP", methValue = "B",
    autoimpute = TRUE, filterDetP = TRUE, ProbeCutoff = 0, SampleCutoff =
    0.1, detPcut = 0.01, filterBeads = TRUE, beadCutoff = 0.05, filterNoCG
    = TRUE, filterSNPs = TRUE, population = NULL, filterMultiHit = TRUE,
    filterXY = TRUE, force = FALSE, arraytype = "450K", verbose = TRUE)
}
\arguments{
  \item{directory}{directory containing the IDAT files (and, for method "minfi", a sample sheet).}
  \item{method}{"ChAMP" (annotation based, default) or "minfi".}
  \item{methValue}{"B" (beta) or "M" (M-value).}
  \item{autoimpute}{impute remaining missing values.}
  \item{filterDetP}{Filter probes and samples on the detection p-value.}
  \item{ProbeCutoff}{Fraction of samples in which a probe may fail before it is removed (0 = any failure removes it).}
  \item{SampleCutoff}{Fraction of failed probes above which a sample is removed.}
  \item{detPcut}{Detection p-value threshold.}
  \item{filterBeads}{Filter probes with bead counts below three.}
  \item{beadCutoff}{Fraction of samples with bead count < 3 above which a probe is removed.}
  \item{filterNoCG}{Remove non-CpG probes.}
  \item{filterSNPs}{Remove probes overlapping known SNPs.}
  \item{population}{SNP mask population.}
  \item{filterMultiHit}{Remove probes that map to multiple genomic locations.}
  \item{filterXY}{Remove probes on the X and Y chromosomes.}
  \item{force}{passed to minfi::read.metharray.exp().}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2" ("EPIC" = "EPICv2").}
  \item{verbose}{print progress.}
}
\description{Load methylation array data}
\value{list with beta (or M), pd, intensity, detP and, for method "minfi", mset and rgSet.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
