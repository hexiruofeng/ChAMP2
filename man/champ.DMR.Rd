\name{champ.DMR}
\alias{champ.DMR}
\title{Differentially methylated regions}
\usage{
champ.DMR(beta = NULL, pheno = NULL, compare.group = NULL, arraytype =
    "450K", method = "Bumphunter", minProbes = 7L, adjPvalDmr = 0.05, cores
    = 3L, maxGap = 300, cutoff = NULL, pickCutoff = TRUE, smooth = TRUE,
    smoothFunction = NULL, useWeights = FALSE, permutations = NULL, B =
    250L, nullMethod = "bootstrap", meanLassoRadius = 375, minDmrSep =
    1000, minDmrSize = 50, adjPvalProbe = 0.05, Rplot = FALSE, PDFplot =
    FALSE, resultsDir = "./CHAMP_ProbeLasso/", rmSNPCH = TRUE, fdr = 0.05,
    dist = 2, mafcut = 0.05, lambda = 1000, C = 2, dmrcate.backend =
    c("auto", "package", "internal"), verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix.}
  \item{pheno}{phenotype vector with exactly two groups (ProbeLasso and DMRcate) or two or more groups (Bumphunter).}
  \item{compare.group}{the two groups to compare when \\code\{pheno\} has more than two levels.}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{method}{"Bumphunter", "ProbeLasso" or "DMRcate".}
  \item{minProbes}{minimum number of probes per region.}
  \item{adjPvalDmr}{adjusted p-value cutoff for reported regions.}
  \item{cores}{number of workers.}
  \item{maxGap}{Maximum distance (bp) between probes in the same cluster.}
  \item{cutoff}{bumphunter cutoff; NULL lets bumphunter pick one.}
  \item{pickCutoff}{Let bumphunter choose the cutoff.}
  \item{smooth}{Smooth the test statistics along the genome.}
  \item{smoothFunction}{Smoothing function passed to bumphunter.}
  \item{useWeights}{Use weights in bumphunter.}
  \item{permutations}{Pre-computed permutations for bumphunter.}
  \item{B}{Number of bumphunter bootstrap permutations.}
  \item{nullMethod}{Null distribution method for bumphunter.}
  \item{meanLassoRadius}{Target mean lasso radius (bp) for ProbeLasso.}
  \item{minDmrSep}{Minimum gap (bp) between merged ProbeLasso regions.}
  \item{minDmrSize}{Minimum DMR width (bp) for ProbeLasso.}
  \item{adjPvalProbe}{Probe-level adjusted p-value cutoff for ProbeLasso.}
  \item{Rplot}{Draw plots on the current device.}
  \item{PDFplot}{Write plots to PDF files in resultsDir.}
  \item{resultsDir}{Directory where result files are written.}
  \item{rmSNPCH}{Remove SNP and cross-hybridising probes before DMRcate.}
  \item{fdr}{Probe-level FDR used by DMRcate.}
  \item{dist}{Distance (bp) used when removing SNP probes for DMRcate.}
  \item{mafcut}{Minor allele frequency cutoff for SNP probe removal.}
  \item{lambda}{Kernel bandwidth (bp) for DMRcate.}
  \item{C}{Bandwidth scaling for DMRcate (sigma = lambda / C).}
  \item{dmrcate.backend}{"auto", "package" or "internal".}
  \item{verbose}{print progress.}
}
\description{Differentially methylated regions}
\value{named list with the DMR table (and, for ProbeLasso, the probe table).}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
