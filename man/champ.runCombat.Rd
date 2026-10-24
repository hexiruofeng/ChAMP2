\name{champ.runCombat}
\alias{champ.runCombat}
\title{Remove batch effects with ComBat}
\usage{
champ.runCombat(beta = NULL, pd = NULL, variablename = "Sample_Group",
    batchname = c("Slide"), logitTrans = TRUE, verbose = TRUE)
}
\arguments{
  \item{beta}{beta value matrix.}
  \item{pd}{sample sheet; must contain \\code\{variablename\} and every entry of \\code\{batchname\}.}
  \item{variablename}{name of the phenotype column that must be protected.}
  \item{batchname}{one or more batch columns to correct.}
  \item{logitTrans}{apply ComBat on the logit scale (recommended for beta values).}
  \item{verbose}{print progress.}
}
\description{Corrects one or more batch variables sequentially, keeping the phenotype of interest (\\code\{variablename\}) in the model so that biological variation is preserved.}
\value{corrected beta value matrix.}
\seealso{\code{\link{champ.process}}, \code{\link{champ_check_env}}}
\keyword{methods}
