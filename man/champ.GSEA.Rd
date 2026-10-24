\name{champ.GSEA}
\alias{champ.GSEA}
\title{Gene set enrichment analysis}
\usage{
champ.GSEA(beta = NULL, DMP = NULL, DMR = NULL, CpGlist = NULL, Genelist =
    NULL, pheno = NULL, method = "fisher", arraytype = "450K", Rplot =
    FALSE, adjPval = 0.05, cores = 1L, verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix (used for annotation and, for "ebayes", the test itself).}
  \item{DMP}{a champ.DMP() result (single data.frame) or NULL.}
  \item{DMR}{a champ.DMR() result or NULL.}
  \item{CpGlist}{character vector or named list of CpG identifiers.}
  \item{Genelist}{character vector or named list of gene symbols.}
  \item{pheno}{phenotype vector (required for method "ebayes").}
  \item{method}{"fisher", "goseq", "gometh" or "ebayes".}
  \item{arraytype}{"450K", "EPICv1" or "EPICv2".}
  \item{Rplot}{draw diagnostic plots.}
  \item{adjPval}{adjusted p-value cutoff.}
  \item{cores}{number of workers for the "ebayes" method.}
  \item{verbose}{print progress.}
}
\description{Four methods are supported: \\describe\{ \\item\{fisher\}\{Fisher's exact test on MSigDB-like gene sets, weighted by the number of CpGs per gene (fast, no extra packages).\} \\item\{goseq\}\{goseq, which models the CpG-count bias explicitly.\} \\item\{gometh\}\{missMethyl::gometh for GO enrichment of CpG lists.\} \\item\{ebayes\}\{champ.ebGSEA(), a global-test based enrichment analysis.\} \}}
\value{named list of enrichment tables.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
