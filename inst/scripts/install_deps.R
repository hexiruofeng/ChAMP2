#!/usr/bin/env Rscript
# =============================================================================
#  ChAMP2 dependency installer
#
#  Installs every package ChAMP can use, on any R from 4.2.3 to 4.6.2 and the
#  matching Bioconductor release (3.16 - 3.24).
#
#  Usage
#  -----
#    Rscript install_deps.R                 # required + optional packages
#    Rscript install_deps.R --required      # required packages only
#    Rscript install_deps.R --check         # only report what is missing
#    Rscript install_deps.R --source        # build from source (needs a toolchain)
#    Rscript install_deps.R --lib=/path     # install into a private library
#
#  The script is intentionally forgiving: ChAMP bundles replacements for the
#  algorithms of the packages that disappeared from CRAN/Bioconductor, so a
#  failure to install one optional package never blocks the core pipeline.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(name) any(grepl(paste0("^--?", name, "(=|$)"), args))
get_value <- function(name, default = NULL) {
    hit <- grep(paste0("^--?", name, "="), args, value = TRUE)
    if (!length(hit)) return(default)
    sub(paste0("^--?", name, "="), "", hit[1])
}

required_only <- get_flag("required") || get_flag("no-suggest")
check_only    <- get_flag("check")
use_source    <- get_flag("source")
lib           <- get_value("lib", NULL)
quiet         <- get_flag("quiet")

msg <- function(...) if (!quiet) cat(sprintf(...), "\n", sep = "")

msg("ChAMP2 dependency installer")
msg("  R version   : %s", R.version.string)
msg("  platform    : %s", R.version$platform)

if (!is.null(lib)) {
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(lib, .libPaths()))
}
msg("  library     : %s", paste(.libPaths(), collapse = ", "))

## ---------------------------------------------------------------------------
## Package manifest
## ---------------------------------------------------------------------------
required <- c(
    "ChAMPdata", "limma", "matrixStats", "illuminaio", "DNAcopy",
    "preprocessCore", "bumphunter", "sva", "globaltest", "GenomicRanges",
    "IRanges", "S4Vectors", "GenomeInfoDb", "quadprog")
optional <- c(
    "minfi", "DMRcate", "impute", "RPMM", "kpmt", "goseq", "geneLenDataBase",
    "missMethyl", "dendextend", "shiny", "shinythemes", "plotly", "DT",
    "prettydoc", "rmarkdown", "knitr", "testthat", "org.Hs.eg.db")

want <- if (required_only) required else c(required, optional)
have <- function(p) nzchar(system.file(package = p))

missing <- want[!vapply(want, have, logical(1))]
msg("  required    : %d/%d installed", sum(vapply(required, have, logical(1))), length(required))
msg("  optional    : %d/%d installed", sum(vapply(optional, have, logical(1))), length(optional))

if (check_only) {
    if (length(missing)) {
        msg("\nMissing packages: %s", paste(missing, collapse = ", "))
        quit(status = if (any(missing %in% required)) 1 else 0)
    }
    msg("\nNothing to install.")
    quit(status = 0)
}
if (!length(missing)) {
    msg("\nEvery requested package is already installed.")
    quit(status = 0)
}
msg("\nInstalling %d package(s): %s", length(missing), paste(missing, collapse = ", "))

## ---------------------------------------------------------------------------
## BiocManager
## ---------------------------------------------------------------------------
if (!have("BiocManager")) {
    msg("\nInstalling BiocManager ...")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (have("BiocManager")) {
    msg("  Bioconductor version: %s", as.character(BiocManager::version()))
}

## A CRAN mirror that also serves extra archived packages.
repos <- c(CRAN = Sys.getenv("CHAMP_CRAN_MIRROR", "https://cloud.r-project.org"))
if (have("BiocManager")) {
    bioc_repos <- tryCatch(BiocManager::repositories(), error = function(e) NULL)
    if (!is.null(bioc_repos)) repos <- c(repos, bioc_repos)
}

type <- if (use_source) "source" else
    if (.Platform$OS.type == "windows") "binary" else
        if (identical(Sys.info()[["sysname"]], "Darwin")) "binary" else "source"

## ---------------------------------------------------------------------------
## Which packages live where
## ---------------------------------------------------------------------------
bioc_known <- c(
    "ChAMPdata", "limma", "illuminaio", "DNAcopy", "preprocessCore",
    "bumphunter", "sva", "globaltest", "GenomicRanges", "IRanges",
    "S4Vectors", "GenomeInfoDb", "minfi", "DMRcate", "impute", "wateRmelon",
    "marray", "goseq", "geneLenDataBase", "missMethyl", "isva", "qvalue",
    "org.Hs.eg.db")

bioc_pkgs <- intersect(missing, bioc_known)
cran_pkgs <- setdiff(missing, bioc_known)

failed <- character(0)

install_bioc <- function(pkgs) {
    if (!length(pkgs)) return(invisible(NULL))
    msg("\n[Bioconductor] %s", paste(pkgs, collapse = ", "))
    if (have("BiocManager")) {
        tryCatch(BiocManager::install(pkgs, update = FALSE, ask = FALSE,
                                      checkBuilt = FALSE, type = type),
                 error = function(e) msg("  BiocManager failed: %s", conditionMessage(e)))
    }
    still <- pkgs[!vapply(pkgs, have, logical(1))]
    if (length(still) && !have("BiocManager")) {
        tryCatch(install.packages(still, repos = repos, type = type),
                 error = function(e) msg("  install.packages failed: %s", conditionMessage(e)))
    }
    invisible(NULL)
}

install_cran <- function(pkgs) {
    if (!length(pkgs)) return(invisible(NULL))
    msg("\n[CRAN] %s", paste(pkgs, collapse = ", "))
    tryCatch(install.packages(pkgs, repos = repos, type = type,
                              dependencies = TRUE, quiet = TRUE),
             error = function(e) msg("  install.packages failed: %s", conditionMessage(e)))
    invisible(NULL)
}

install_bioc(bioc_pkgs)
install_cran(cran_pkgs)

## ---------------------------------------------------------------------------
## Fallbacks for packages that are not (or no longer) on CRAN / Bioconductor
## ---------------------------------------------------------------------------
have_all <- function(p) vapply(p, have, logical(1))

if ("ChAMPdata" %in% missing && !have("ChAMPdata")) {
    msg("\nChAMPdata is not available from the configured repositories; trying GitHub ...")
    url <- "https://github.com/YuanTian1991/ChAMPdata/archive/refs/heads/master.tar.gz"
    tryCatch({
        tmp <- file.path(tempdir(), "ChAMPdata.tar.gz")
        utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
        install.packages(tmp, repos = NULL, type = "source")
    }, error = function(e) msg("  GitHub install failed: %s", conditionMessage(e)))
}

if ("kpmt" %in% missing && !have("kpmt")) {
    msg("\nNote: kpmt was archived on CRAN (2025-04-27). ChAMP bundles champ_kpmt(),")
    msg("      so the KPMT column of champ.ebGSEA() works without it.")
}

if ("RPMM" %in% missing && !have("RPMM")) {
    msg("\nNote: ChAMP bundles champ_blc(), so BMIQ works without RPMM.")
}
if ("impute" %in% missing && !have("impute")) {
    msg("Note: ChAMP bundles champ_knn_impute(), so imputation works without impute.")
}
if ("DMRcate" %in% missing && !have("DMRcate")) {
    msg("Note: ChAMP contains a DMRcate-compatible implementation, so method = 'DMRcate' works.")
}

## ---------------------------------------------------------------------------
## Report
## ---------------------------------------------------------------------------
still_missing <- want[!vapply(want, have, logical(1))]
msg("\n--------------------------------------------")
msg("Installed : %d/%d", sum(vapply(want, have, logical(1))), length(want))
if (length(still_missing)) {
    msg("Missing   : %s", paste(still_missing, collapse = ", "))
    msg("Required packages missing: %s",
        paste(intersect(still_missing, required), collapse = ", "))
    if (length(setdiff(still_missing, required))) {
        msg("ChAMP still works; the features that need the missing optional packages")
        msg("will report exactly what to install.")
    }
} else {
    msg("All requested packages are installed.")
}
msg("--------------------------------------------")

status <- if (all(vapply(required, have, logical(1)))) 0 else 1
quit(status = status)