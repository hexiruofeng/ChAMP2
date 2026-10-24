# =============================================================================
#  Dependency management
#
#  ChAMP 2.32 keeps a single, machine readable dependency manifest so that
#  champ_install_deps() and champ_check_env() can never drift away from
#  DESCRIPTION.  Optional packages are *packages*, not hard requirements:
#  every feature that needs one fails with an actionable message instead of a
#  cryptic "there is no package called ..." at load time.
# =============================================================================

#' ChAMP dependency manifest
#'
#' @param which "Imports", "Suggests" or "all".
#' @return data.frame with package, role, source, min_version and feature.
#' @export
champ_dependencies <- function(which = c("all", "Imports", "Suggests")) {
    which <- match.arg(which)
    man <- .champ_dependency_manifest()
    if (!identical(which, "all")) man <- man[man$role == which, , drop = FALSE]
    rownames(man) <- NULL
    man
}

#' Report the status of every ChAMP dependency
#'
#' Prints (and invisibly returns) a table with the installed version of every
#' package ChAMP can use, whether it is required, and which feature it enables.
#'
#' @param quiet do not print the report.
#' @return invisibly, a data.frame.
#' @export
champ_check_env <- function(quiet = FALSE) {
    tab <- .champ_dependency_table()
    tab$ok <- ifelse(tab$installed, "ok", ifelse(tab$role == "Imports", "MISSING", "optional"))
    if (!isTRUE(quiet)) {
        cat("ChAMP dependency report\n")
        cat("  R          : ", R.version.string, "\n", sep = "")
        cat("  Bioc       : ", .champ_bioc_version(), "\n", sep = "")
        cat("  library    : ", paste(.libPaths(), collapse = ", "), "\n\n", sep = "")
        print(tab[, c("package", "role", "installed", "version", "ok")], row.names = FALSE)
        miss <- tab$package[!tab$installed & tab$role == "Imports"]
        if (length(miss)) {
            cat("\nRequired packages missing: ", paste(miss, collapse = ", "), "\n", sep = "")
            cat("Run ChAMP::champ_install_deps() to install them.\n")
        } else {
            cat("\nAll required packages are installed.\n")
        }
        opt <- tab$package[!tab$installed & tab$role == "Suggests"]
        if (length(opt)) {
            cat("Optional packages not installed: ", paste(opt, collapse = ", "), "\n", sep = "")
        }
    }
    invisible(tab)
}

#' Current Bioconductor version
#' @keywords internal
.champ_bioc_version <- function() {
    if (.champ_has("BiocManager")) {
        return(.champ_try(as.character(BiocManager::version()), default = "unknown"))
    }
    "unknown (BiocManager not installed)"
}

#' Repositories used by champ_install_deps()
#' @keywords internal
.champ_repos <- function() {
    repos <- c(CRAN = "https://cloud.r-project.org")
    if (.champ_has("BiocManager")) {
        bioc <- .champ_try(BiocManager::repositories(), default = NULL)
        if (!is.null(bioc)) repos <- c(repos, bioc)
    }
    repos
}

#' Install the packages ChAMP needs
#'
#' Installs every missing required package (Bioconductor first, then CRAN) and,
#' unless \code{optional = FALSE}, the suggested packages that enable the
#' optional features (minfi loading, DMRcate backend, RPMM/impute/kpmt legacy
#' backends, the Shiny GUIs, ...).
#'
#' The function is deliberately forgiving: a package that cannot be installed
#' (for example because the current Bioconductor release dropped it) produces a
#' warning and the installation continues, because ChAMP bundles a fallback for
#' each of those algorithms.
#'
#' @param optional also install the Suggests packages.
#' @param binary use binary packages when available (much faster on Windows and
#'   macOS).
#' @param upgrade "never" (default), "default" or "always".
#' @param force reinstall packages that are already present.
#' @param verbose print progress.
#'
#' @return invisibly, a list with installed, skipped and failed packages.
#' @export
champ_install_deps <- function(optional = TRUE, binary = TRUE, upgrade = "never",
                               force = FALSE, verbose = TRUE) {
    if (!.champ_has("BiocManager")) {
        .champ_msg("Installing BiocManager...")
        utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
        if (!.champ_has("BiocManager")) {
            .champ_stop("Could not install BiocManager; check your network/proxy settings.")
        }
    }
    repos <- .champ_repos()
    type <- if (isTRUE(binary) && .Platform$OS.type == "windows") "binary" else
        if (isTRUE(binary) && identical(Sys.info()[["sysname"]], "Darwin")) "binary" else "source"

    man <- .champ_dependency_manifest()
    if (!isTRUE(optional)) man <- man[man$role == "Imports", , drop = FALSE]
    need <- if (isTRUE(force)) man$package else
        man$package[!vapply(man$package, .champ_has, logical(1L))]
    if (!length(need)) {
        .champ_msg("All requested ChAMP dependencies are already installed.")
        return(invisible(list(installed = character(0), skipped = man$package,
                              failed = character(0))))
    }

    bioc <- intersect(need, man$package[man$source == "Bioconductor"])
    cran <- setdiff(need, bioc)
    installed <- character(0)
    failed <- character(0)

    if (length(bioc)) {
        .champ_msg("Installing ", length(bioc), " Bioconductor package(s): ",
                   paste(bioc, collapse = ", "))
        ok <- .champ_try({
            BiocManager::install(bioc, update = FALSE, ask = FALSE,
                                 checkBuilt = FALSE, type = type)
            TRUE
        }, default = FALSE)
        if (!ok) .champ_warn("BiocManager::install() reported a problem; continuing.")
        installed <- c(installed, bioc[vapply(bioc, .champ_has, logical(1L))])
    }

    if (length(cran)) {
        .champ_msg("Installing ", length(cran), " CRAN package(s): ",
                   paste(cran, collapse = ", "))
        .champ_try(utils::install.packages(cran, repos = repos, type = type,
                                           dependencies = TRUE, quiet = !verbose),
                   default = NULL)
        installed <- c(installed, cran[vapply(cran, .champ_has, logical(1L))])
    }

    ## kpmt was archived on CRAN (2025-04-27); install the last source release
    ## when it is requested and the fallback is not enough for the user.
    if ("kpmt" %in% need && !.champ_has("kpmt")) {
        .champ_msg("kpmt is no longer on CRAN; ChAMP ships champ.kpmt() as a replacement.")
    }

    failed <- setdiff(need, installed)
    if (length(failed)) {
        .champ_warn("These packages could not be installed: ", paste(failed, collapse = ", "),
                    "\n  ChAMP remains usable for everything that does not need them; ",
                    "run ChAMP::champ_check_env() for details.")
    }
    .champ_msg("Done. ", length(installed), " package(s) installed, ",
               length(failed), " failed.")
    invisible(list(installed = installed,
                   skipped = setdiff(man$package, need),
                   failed = failed))
}
