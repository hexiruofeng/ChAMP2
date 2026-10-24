\name{champ.filter}
\alias{champ.filter}
\title{Filter a methylation matrix}
\usage{
champ.filter(beta = NULL, M = NULL, pd = NULL, intensity = NULL, Meth =
    NULL, UnMeth = NULL, detP = NULL, beadcount = NULL, autoimpute = TRUE,
    filterDetP = TRUE, ProbeCutoff = 0, SampleCutoff = 0.1, detPcut = 0.01,
    filterBeads = TRUE, beadCutoff = 0.05, filterNoCG = TRUE, filterSNPs =
    TRUE, population = NULL, filterMultiHit = TRUE, filterXY = TRUE,
    fixOutlier = TRUE, arraytype = "450K", verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix (probes x samples); can be NULL when M is given.}
  \item{M}{M value matrix; can be NULL when beta is given.}
  \item{pd}{sample sheet.}
  \item{intensity}{A matrix of total intensities (meth + unmeth), probes x samples.}
  \item{Meth}{Optional matrix of methylated intensities, filtered alongside the data.}
  \item{UnMeth}{Optional matrix of unmethylated intensities, filtered alongside the data.}
  \item{detP}{detection p-value matrix.}
  \item{beadcount}{bead count matrix.}
  \item{autoimpute}{impute remaining missing values (KNN, k = 5).}
  \item{filterDetP}{Filter probes and samples on the detection p-value.}
  \item{ProbeCutoff}{fraction of samples that may fail for a probe (0 removes any probe failing in one sample).}
  \item{SampleCutoff}{fraction of failed probes above which a sample is removed.}
  \item{detPcut}{detection p-value threshold.}
  \item{filterBeads}{Filter probes with bead counts below three.}
  \item{beadCutoff}{fraction of samples with bead count < 3 above which a probe is removed.}
  \item{filterNoCG}{Remove non-CpG probes.}
  \item{filterSNPs}{Remove probes overlapping known SNPs.}
  \item{population}{SNP mask population (e.g. "EUR", "AFR"); NULL uses the general mask.}
  \item{filterMultiHit}{Remove probes that map to multiple genomic locations.}
  \item{filterXY}{Remove probes on the X and Y chromosomes.}
  \item{fixOutlier}{replace exact 0/1 beta values.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{verbose}{print progress.}
}
\description{Applies ChAMP's standard filtering steps to a beta (or M) value matrix: detection p-value filtering, optional KNN imputation, bead count filtering, non-CpG removal, SNP masking, multi-hit removal and sex chromosome removal.}
\value{list of the filtered matrices plus pd.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
