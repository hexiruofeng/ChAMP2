\name{champnews}
\title{ChAMP News}
\encoding{UTF-8}

\section{Version 2.32.0}{
  \itemize{
    \item One source tree now supports R 4.2.3 to 4.6.2 (Bioconductor 3.16 to 3.24).
    \item Hard dependencies reduced to 14 packages; minfi, DMRcate, DT, RPMM,
          impute, marray, wateRmelon and Illumina450ProbeVariants.db are optional.
    \item Bundled algorithms: champ_knn_impute() (impute), champ_blc() (RPMM),
          champ_kpmt() (kpmt) and champ_dmrcate() (DMRcate).
    \item New champ_check_env(), champ_install_deps(), champ_dependencies() and
          inst/scripts/install_deps.R.
    \item Annotation loading centralised in .champ_probe_features(); EPICv1 and
          EPICv2 work end to end, ProbeLasso no longer needs illumina450Gr.
    \item Fixed champ.import(), champ.filter(), champ.process(), champ.BMIQ(),
          champ.norm(), champ.SVD(), champ.DMP(), champ.DMR(), champ.CNA(),
          champ.refbase() and champ.runCombat() bugs listed in NEWS.md.
    \item champ.QC() no longer needs minfi; parallel execution uses base parallel.
  }
}

\section{Version 2.29.1}{
  \itemize{
    \item EPICv2 support; \code{arraytype = "EPIC"} now means EPICv2.
  }
}

\section{Version 2.18.2}{
  \itemize{
    \item Added parallel running for ebGSEA function and the enriched gene list.
  }
}

\section{Version 2.17.3}{
  \itemize{
    \item Fixed a bug in runCombat when some factors are numeric.
  }
}
