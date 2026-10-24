#!/usr/bin/env Rscript
# =============================================================================
#  bundle_deps.R -- create an offline dependency repository for ChAMP2
#
#  Downloads every package ChAMP needs (and, optionally, the suggested ones)
#  as binaries for the running R/Bioconductor version into a local repository
#  that can be copied to a machine without internet access:
#
#      Rscript tools/bundle_deps.R --dest=D:/champ2-deps
#      # on the offline machine:
#      Rscript tools/bundle_deps.R --install --dest=D:/champ2-deps
#
#  The repository layout is the standard CRAN/Bioconductor one
#  (src/contrib or bin/windows/contrib/<R.x>), so install.packages(repos =
#  "file:///D:/champ2-deps", type = "binary") also works.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
flag <- function(nm) any(grepl(paste0("^--?", nm, "(=|$)"), args))
value <- function(nm, default = NULL) {
    hit <- grep(paste0("^--?", nm, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^--?", nm, "="), "", hit[1])
}

dest <- value("dest", file.path(getwd(), "champ2-deps"))
install_mode <- flag("install")

required <- c("ChAMPdata", "limma", "matrixStats", "illuminaio", "DNAcopy",
              "preprocessCore", "bumphunter", "sva", "globaltest", "GenomicRanges",
              "IRanges", "S4Vectors", "GenomeInfoDb", "quadprog")
optional <- c("minfi", "DMRcate", "impute", "RPMM", "goseq", "geneLenDataBase",
              "missMethyl", "dendextend", "ggplot2", "shiny", "shinythemes",
              "plotly", "DT", "knitr", "rmarkdown", "prettydoc", "testthat")
want <- if (flag("required")) required else c(required, optional)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

subdir <- if (.Platform$OS.type == "windows") {
    file.path("bin", "windows", "contrib", paste0(R.version$major, ".", R.version$minor))
} else {
    file.path("src", "contrib")
}
repo <- file.path(dest, subdir)
dir.create(repo, recursive = TRUE, showWarnings = FALSE)

if (install_mode) {
    message("Installing from the local repository ", dest)
    install.packages(want, repos = c(CRAN = paste0("file:///", normalizePath(dest))),
                     type = if (.Platform$OS.type == "windows") "binary" else "source")
    quit(status = 0)
}

message("Downloading ", length(want), " packages into ", repo)
utils::download.packages(want, destdir = repo,
                         repos = BiocManager::repositories(),
                         type = if (.Platform$OS.type == "windows") "binary" else "source")
message("Done. Copy ", dest, " to the offline machine and run:")
message("  Rscript tools/bundle_deps.R --install --dest=", dest)