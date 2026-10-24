\name{champ_dmrcate}
\alias{champ_dmrcate}
\title{DMRcate-compatible DMR calling}
\usage{
champ_dmrcate(beta, pheno, arraytype = "450K", fdr = 0.05, min.cpgs = 7L,
    lambda = 1000, C = 2, pcutoff = "fdr", betacutoff = NULL, verbose =
    TRUE)
}
\arguments{
  \item{beta}{beta value matrix (probes x samples).}
  \item{pheno}{phenotype vector with exactly two groups (or a numeric covariate).}
  \item{arraytype}{one of "450K", "EPICv1", "EPICv2".}
  \item{fdr}{probe-level FDR used to define significant probes.}
  \item{min.cpgs}{minimum number of probes per DMR.}
  \item{lambda}{kernel bandwidth in base pairs.}
  \item{C}{scaling of the bandwidth (\\code\{sigma = lambda / C\}).}
  \item{pcutoff}{"fdr" or a numeric probe-level p-value cutoff.}
  \item{betacutoff}{optional absolute mean beta difference filter.}
  \item{verbose}{print progress.}
}
\description{A self-contained implementation of the DMRcate algorithm (Peters et al., Epigenetics & Chromatin 2015): probes are annotated with moderated t-statistics, a Gaussian kernel with bandwidth \\code\{lambda / C\} is used to compute local chi-square statistics, and regions of consecutive significant probes are summarised with Stouffer, harmonic-mean and Fisher combination of p-values. When the \\code\{DMRcate\} package is installed, \\code\{champ.DMR()\} uses it by default; this function is the fallback, and produces a data frame with the same columns as \\code\{as.data.frame(DMRcate::extractRanges(...))\}.}
\value{data.frame with columns seqnames, start, end, width, strand, no.cpgs, min_smoothed_fdr, Stouffer, HMFDR, Fisher, maxdiff, meandiff, overlapping.genes.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
