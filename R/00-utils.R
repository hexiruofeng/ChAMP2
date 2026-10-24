# =============================================================================
#  ChAMP 2.32 -- internal utilities and compatibility layer
#
#  This file is the foundation of the R 4.2.3 -- 4.6.2 compatibility work.  It
#  provides:
#    * small helpers used by every exported function;
#    * a single, cached access point for the ChAMPdata annotation objects (the
#      old code called utils::data() in about twenty places with hard-coded
#      object names, several of which no longer exist in recent ChAMPdata);
#    * array-type normalisation (450K / EPICv1 / EPICv2 / 27K / Mouse);
#    * pure-R replacements for beta/m-value transformations so that functions
#      no longer need minfi just to call logit2()/ilogit2().
#
#  Nothing in this file is exported.
# =============================================================================

# ---------------------------------------------------------------------------
# package state
# ---------------------------------------------------------------------------

.champ_state <- new.env(parent = emptyenv())
.champ_state$cache  <- new.env(parent = emptyenv())
.champ_state$warned <- new.env(parent = emptyenv())

.champ_version <- function() {
    as.character(utils::packageVersion("ChAMP"))
}

.champ_or <- function(x, y) if (is.null(x)) y else x

.champ_msg <- function(...) {
    if (isTRUE(getOption("ChAMP.quiet", FALSE))) return(invisible(NULL))
    message(...)
    invisible(NULL)
}

.champ_stop <- function(...) {
    stop(paste0(...), call. = FALSE)
}

.champ_warn <- function(...) {
    warning(paste0(...), call. = FALSE, immediate. = TRUE)
    invisible(NULL)
}

.champ_warn_once <- function(key, ...) {
    if (!is.null(.champ_state$warned[[key]])) return(invisible(NULL))
    .champ_state$warned[[key]] <- TRUE
    .champ_warn(...)
}

# ---------------------------------------------------------------------------
# optional package handling
# ---------------------------------------------------------------------------

.champ_has <- function(pkg) {
    isTRUE(requireNamespace(pkg, quietly = TRUE))
}

.champ_require <- function(pkg, feature = NULL, hint = TRUE) {
    if (.champ_has(pkg)) return(invisible(TRUE))
    what <- if (is.null(feature)) "" else paste0(" for ", feature)
    msg <- sprintf("Package '%s' is required%s but is not installed.", pkg, what)
    if (isTRUE(hint)) {
        msg <- paste0(
            msg, "\n",
            "  * Install all ChAMP dependencies at once with:\n",
            "      ChAMP::champ_install_deps()\n",
            "  * or install just this package with:\n",
            "      BiocManager::install(\"", pkg, "\")\n",
            "  * see ChAMP::champ_check_env() for a full dependency report."
        )
    }
    .champ_stop(msg)
}

.champ_try <- function(expr, default = NULL) {
    tryCatch(expr, error = function(e) default)
}

# ---------------------------------------------------------------------------
# array type handling
# ---------------------------------------------------------------------------

.champ_array_types <- c("450K", "EPIC", "EPICv1", "EPICv2", "27K", "Mouse")

#' Normalise a user supplied array type.
#'
#' Since ChAMP 2.29 "EPIC" means EPIC v2.  We keep that behaviour (it is what
#' upstream documents) but accept every spelling seen in the wild and make the
#' mapping explicit instead of silently guessing.
.champ_arraytype <- function(arraytype = "450K") {
    if (is.null(arraytype) || !length(arraytype) || is.na(arraytype[1L])) {
        .champ_stop("'arraytype' must be one of: ", paste(.champ_array_types, collapse = ", "))
    }
    a <- toupper(trimws(as.character(arraytype)[1L]))
    a <- gsub("[ _-]", "", a)
    out <- switch(a,
        "450K"   = "450K",
        "HM450"  = "450K",
        "450"    = "450K",
        "ILLUMINAHUMANMETHYLATION450K" = "450K",
        "EPIC"   = "EPICv2",
        "EPICV2" = "EPICv2",
        "EPIC2"  = "EPICv2",
        "V2"     = "EPICv2",
        "EPICV1" = "EPICv1",
        "EPIC1"  = "EPICv1",
        "V1"     = "EPICv1",
        "27K"    = "27K",
        "HM27"   = "27K",
        "MOUSE"  = "Mouse",
        NA_character_)
    if (is.na(out)) {
        .champ_stop("Unsupported 'arraytype': '", arraytype, "'. Supported values are: ",
                    paste(.champ_array_types, collapse = ", "), ".")
    }
    if (identical(a, "EPIC")) {
        .champ_warn_once("epic_is_v2",
            "arraytype = 'EPIC' is interpreted as EPIC version 2 (the ChAMP >= 2.29 default). ",
            "Use arraytype = 'EPICv1' for the original (ilm10b4.hg19) EPIC array.")
    }
    out
}

#' Genome build implied by an array type.
.champ_genome <- function(arraytype = "450K") {
    at <- .champ_arraytype(arraytype)
    switch(at,
        "450K"   = "hg19",
        "EPICv1" = "hg19",
        "EPICv2" = "hg38",
        "27K"    = "hg19",
        "Mouse"  = "mm10")
}

# ---------------------------------------------------------------------------
# ChAMPdata / annotation access
# ---------------------------------------------------------------------------

#' Load a data object from a package, robust to name != file name.
#'
#' utils::data() only works when the requested name matches an entry of the
#' package's data/datalist; several ChAMPdata objects (most notably bloodCtl,
#' stored in champBloodCtls.rda) are not reachable that way, and some objects
#' were renamed between ChAMPdata releases.  This helper tries, in order: the
#' calling environment, the package namespace, a direct load() of a matching
#' .rda file, and finally utils::data().
.champ_data_object <- function(name, package = "ChAMPdata", required = TRUE,
                               envir = parent.frame(), files = NULL) {
    if (exists(name, envir = envir, inherits = FALSE)) {
        return(get(name, envir = envir))
    }
    if (!.champ_has(package)) {
        if (isTRUE(required)) {
            .champ_stop("Package '", package, "' is required to load '", name, "' but is not installed.\n",
                        "  Install it with ChAMP::champ_install_deps() or BiocManager::install(\"", package, "\").")
        }
        return(NULL)
    }
    ns <- asNamespace(package)
    if (exists(name, envir = ns, inherits = FALSE)) {
        obj <- get(name, envir = ns)
        assign(name, obj, envir = envir)
        return(obj)
    }

    data_dir <- system.file("data", package = package)
    if (nzchar(data_dir)) {
        cand <- unique(c(files,
                         paste0(name, c(".rda", ".RData", ".rdata")),
                         paste0(gsub("\\.", "_", name), c(".rda", ".RData"))))
        for (f in cand) {
            p <- file.path(data_dir, f)
            if (file.exists(p)) {
                e <- new.env(parent = emptyenv())
                ok <- tryCatch({ load(p, envir = e); TRUE }, error = function(err) FALSE)
                if (ok && length(ls(e, all.names = TRUE))) {
                    nms <- ls(e, all.names = TRUE)
                    for (nm in nms) assign(nm, get(nm, envir = e), envir = envir)
                    if (name %in% nms) return(get(name, envir = envir))
                    return(get(nms[1L], envir = envir))
                }
            }
        }
    }

    e <- new.env(parent = emptyenv())
    ok <- tryCatch({
        utils::data(list = name, package = package, envir = e)
        length(ls(e, all.names = TRUE)) > 0L
    }, error = function(err) FALSE,
       warning = function(w) length(ls(e, all.names = TRUE)) > 0L)
    if (ok) {
        nms <- ls(e, all.names = TRUE)
        for (nm in nms) assign(nm, get(nm, envir = e), envir = envir)
        if (name %in% nms) return(get(name, envir = envir))
        return(get(nms[1L], envir = envir))
    }
    if (isTRUE(required)) {
        .champ_stop("Could not load dataset '", name, "' from package '", package, "'.\n",
                    "  Your ", package, " installation may be too old or corrupted; please reinstall it with\n",
                    "      ChAMP::champ_install_deps()")
    }
    NULL
}

# -- small normalisers for annotation columns --------------------------------

.champ_norm_chr <- function(x) {
    x <- toupper(trimws(as.character(x)))
    x <- sub("^CHR", "", x)
    x[x %in% "23"] <- "X"
    x[x %in% "24"] <- "Y"
    x
}

.champ_norm_cgi <- function(x) {
    y <- tolower(gsub("[ _.-]", "", trimws(as.character(x))))
    y <- sub("^open$", "opensea", y)
    y[grepl("^opensea", y)] <- "opensea"
    y[grepl("^island", y)]  <- "island"
    y[grepl("^shelf", y)]   <- "shelf"
    y[grepl("^shore", y)]   <- "shore"
    y
}

.champ_norm_feature <- function(x) {
    y <- trimws(as.character(x))
    y[tolower(y) %in% c("1stexon", "1st exon", "firstexon")] <- "1stExon"
    y[tolower(y) %in% c("3'utr", "3utr", "3' utr")] <- "3'UTR"
    y[tolower(y) %in% c("5'utr", "5utr", "5' utr")] <- "5'UTR"
    y[tolower(y) == "body"]   <- "Body"
    y[tolower(y) == "igr"]    <- "IGR"
    y[tolower(y) == "tss1500"] <- "TSS1500"
    y[tolower(y) == "tss200"]  <- "TSS200"
    y
}

# -- probe annotation --------------------------------------------------------

.champ_probe_candidates <- function(at) {
    switch(at,
        "450K"   = c("probe.features"),
        "EPICv1" = c("probe.features.epicv1", "probe.features.epic"),
        "EPICv2" = c("probe.features.epicv2", "probe.features.epic"),
        "27K"    = c("probe.features"),
        "Mouse"  = c("probe.features.mouse"),
        character(0))
}

#' Coerce any probe annotation table to the ChAMP layout.
.champ_prepare_probe_features <- function(pf, arraytype = "450K") {
    pf <- as.data.frame(pf, stringsAsFactors = FALSE)
    if (is.null(rownames(pf)) || all(rownames(pf) == as.character(seq_len(nrow(pf))))) {
        idcol <- intersect(c("Name", "probeID", "IlmnID", "ID"), colnames(pf))
        if (length(idcol)) rownames(pf) <- as.character(pf[[idcol[1L]]])
    }
    if (is.null(pf$Name)) pf$Name <- rownames(pf)

    chr_col <- intersect(c("CHR", "chr", "Chromosome", "seqnames", "chr_name"), colnames(pf))
    if (!length(chr_col)) .champ_stop("Probe annotation has no chromosome column.")
    pf$CHR <- .champ_norm_chr(pf[[chr_col[1L]]])

    pos_col <- intersect(c("MAPINFO", "pos", "position", "start", "POS"), colnames(pf))
    if (!length(pos_col)) .champ_stop("Probe annotation has no position column.")
    pf$MAPINFO <- suppressWarnings(as.numeric(pf[[pos_col[1L]]]))

    cgi_col <- intersect(c("cgi", "Relation_to_Island", "relation_to_island", "CGI"), colnames(pf))
    pf$cgi <- if (length(cgi_col)) .champ_norm_cgi(pf[[cgi_col[1L]]]) else NA_character_

    fea_col <- intersect(c("feature", "UCSC_RefGene_Group", "Feature"), colnames(pf))
    pf$feature <- if (length(fea_col)) .champ_norm_feature(pf[[fea_col[1L]]]) else NA_character_

    if (is.null(pf$gene)) {
        gcol <- intersect(c("gene", "UCSC_RefGene_Name", "Gene_Symbol", "symbol"), colnames(pf))
        pf$gene <- if (length(gcol)) as.character(pf[[gcol[1L]]]) else ""
    } else {
        pf$gene <- as.character(pf$gene)
    }
    pf$gene[is.na(pf$gene)] <- ""

    if (is.null(pf$Type)) {
        tcol <- intersect(c("Type", "type", "Infinium_Design_Type"), colnames(pf))
        pf$Type <- if (length(tcol)) as.character(pf[[tcol[1L]]]) else NA_character_
    }
    pf$Type <- as.character(pf$Type)
    t1 <- !is.na(pf$Type) & tolower(pf$Type) %in% c("i", "1", "type i", "infinium i")
    pf$Type[t1] <- "I"
    t2 <- !is.na(pf$Type) & tolower(pf$Type) %in% c("ii", "2", "type ii", "infinium ii")
    pf$Type[t2] <- "II"

    if (is.null(pf$feat.cgi)) {
        pf$feat.cgi <- paste(pf$feature, pf$cgi, sep = " - ")
    }
    pf
}

#' Build a probe annotation table from a minfi annotation package.
.champ_probe_features_minfi <- function(at) {
    ann <- switch(at,
        "450K"   = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        "EPICv1" = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
        "EPICv2" = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
        "27K"    = "IlluminaHumanMethylation27kanno.ilmn12.hg19",
        NULL)
    if (is.null(ann) || !.champ_has(ann) || !.champ_has("minfi")) return(NULL)
    anno <- .champ_try({
        obj <- getExportedValue(ann, ann)
        as.data.frame(minfi::getAnnotation(obj))
    })
    if (is.null(anno) || !nrow(anno)) return(NULL)
    anno
}

#' Cached probe annotation for an array type.
.champ_probe_features <- function(arraytype = "450K") {
    at <- .champ_arraytype(arraytype)
    key <- paste0("pf:", at)
    hit <- .champ_state$cache[[key]]
    if (!is.null(hit)) return(hit)

    pf <- NULL
    src <- NA_character_
    for (nm in .champ_probe_candidates(at)) {
        obj <- .champ_try(.champ_data_object(nm, required = FALSE), default = NULL)
        if (is.data.frame(obj) && nrow(obj) > 0L) {
            pf <- obj
            src <- nm
            break
        }
    }
    if (is.null(pf)) {
        pf <- .champ_probe_features_minfi(at)
        src <- if (is.null(pf)) NA_character_ else "minfi"
    }
    if (is.null(pf)) {
        .champ_stop(
            "No probe annotation available for arraytype = '", at, "'.\n",
            "  ChAMP needs the ChAMPdata package (>= 2.6.0).  Install it with:\n",
            "      ChAMP::champ_install_deps()\n",
            "  or: BiocManager::install(\"ChAMPdata\")")
    }
    pf <- .champ_prepare_probe_features(pf, at)
    attr(pf, "source")    <- src
    attr(pf, "arraytype") <- at
    attr(pf, "genome")    <- .champ_genome(at)
    .champ_state$cache[[key]] <- pf
    pf
}

#' Convenience: annotation restricted to a set of probes, in the requested order.
.champ_probe_features_for <- function(probes, arraytype = "450K") {
    pf <- .champ_probe_features(arraytype)
    probes <- as.character(probes)
    out <- pf[match(probes, rownames(pf)), , drop = FALSE]
    rownames(out) <- probes
    out
}

# -- SNP masks ---------------------------------------------------------------

.champ_probe_mask <- function(arraytype = "450K") {
    at <- .champ_arraytype(arraytype)
    key <- paste0("mask:", at)
    hit <- .champ_state$cache[[key]]
    if (!is.null(hit)) return(hit)

    cand <- switch(at,
        "450K"   = list("Hm450K_Mask", "hm450.manifest.hg19"),
        "EPICv1" = list("EPICv1_Mask", "EPIC.manifest.hg19"),
        "EPICv2" = list("EPICv2_Mask"),
        "27K"    = list(),
        "Mouse"  = list())

    mask <- NULL
    for (nm in cand) {
        obj <- .champ_try(.champ_data_object(nm, required = FALSE), default = NULL)
        if (is.data.frame(obj) && nrow(obj) > 0L) { mask <- obj; break }
    }
    if (is.null(mask)) return(NULL)

    mask <- as.data.frame(mask, stringsAsFactors = FALSE)
    idcol <- intersect(c("probeID", "Name", "probe"), colnames(mask))
    if (length(idcol) && !identical(rownames(mask), as.character(mask[[idcol[1L]]]))) {
        rownames(mask) <- as.character(mask[[idcol[1L]]])
    }
    .champ_state$cache[[key]] <- mask
    mask
}

#' Masked probe identifiers for a population ("general").
.champ_mask_probes <- function(arraytype = "450K", population = NULL) {
    mask <- .champ_probe_mask(arraytype)
    if (is.null(mask)) return(character(0))
    if (!is.null(population) && length(population) == 1L) {
        col <- paste0("MASK_general_", population)
        if (col %in% colnames(mask)) {
            .champ_msg("    Using ", population, " specific SNP mask.")
            return(rownames(mask)[which(as.logical(mask[[col]]))])
        }
        .champ_msg("    Population '", population, "' is not available for this array; using the general mask.")
    }
    if ("MASK_general" %in% colnames(mask)) {
        return(rownames(mask)[which(as.logical(mask[["MASK_general"]]))])
    }
    character(0)
}

# -- other ChAMPdata objects -------------------------------------------------

.champ_multi_hit <- function() {
    hit <- .champ_state$cache[["multi.hit"]]
    if (!is.null(hit)) return(hit)
    mh <- .champ_try(.champ_data_object("multi.hit", required = FALSE), default = NULL)
    if (is.null(mh)) {
        mh <- .champ_try(.champ_data_object("multi.hit", files = "multi.hit.rda",
                                            required = FALSE), default = NULL)
    }
    if (is.null(mh)) return(NULL)
    mh <- as.data.frame(mh, stringsAsFactors = FALSE)
    idcol <- intersect(c("TargetID", "probeID", "Name"), colnames(mh))
    if (!length(idcol)) return(NULL)
    ids <- unique(as.character(mh[[idcol[1L]]]))
    .champ_state$cache[["multi.hit"]] <- ids
    ids
}

#' MSigDB-like pathway list used by champ.GSEA / champ.ebGSEA.
.champ_pathway_list <- function(required = TRUE) {
    hit <- .champ_state$cache[["PathwayList"]]
    if (!is.null(hit)) return(hit)
    pl <- .champ_try(.champ_data_object("PathwayList", required = FALSE), default = NULL)
    if (is.null(pl)) {
        if (isTRUE(required)) {
            .champ_stop("ChAMPdata::PathwayList is required for gene-set enrichment but could not be loaded.\n",
                        "  Reinstall ChAMPdata with ChAMP::champ_install_deps().")
        }
        return(NULL)
    }
    pl <- as.list(pl)
    pl <- pl[!vapply(pl, is.null, logical(1L))]
    .champ_state$cache[["PathwayList"]] <- pl
    pl
}

.champ_control_probes <- function(arraytype = "450K") {
    at <- .champ_arraytype(arraytype)
    key <- paste0("ctl:", at)
    hit <- .champ_state$cache[[key]]
    if (!is.null(hit)) return(hit)
    nm <- if (at %in% c("EPIC", "EPICv1", "EPICv2")) "ControlProbesEPIC" else "ControlProbes450K"
    cp <- .champ_try(.champ_data_object(nm, required = FALSE), default = NULL)
    if (is.null(cp)) return(NULL)
    .champ_state$cache[[key]] <- as.data.frame(cp, stringsAsFactors = FALSE)
    .champ_state$cache[[key]]
}

.champ_blood_controls <- function() {
    hit <- .champ_state$cache[["bloodCtl"]]
    if (!is.null(hit)) return(hit)
    bc <- .champ_try(.champ_data_object("bloodCtl", required = FALSE), default = NULL)
    if (is.null(bc)) {
        bc <- .champ_try(.champ_data_object("champBloodCtls", files = "champBloodCtls.rda",
                                            required = FALSE), default = NULL)
    }
    if (is.null(bc)) return(NULL)
    .champ_state$cache[["bloodCtl"]] <- bc
    bc
}

.champ_cell_type_means <- function(arraytype = "450K") {
    at <- .champ_arraytype(arraytype)
    nm <- if (identical(at, "27K")) "CellTypeMeans27K" else "CellTypeMeans450K"
    .champ_try(.champ_data_object(nm, required = FALSE), default = NULL)
}

.champ_match_gene_name <- function() {
    .champ_try(.champ_data_object("MatchGeneName", required = FALSE), default = NULL)
}

# ---------------------------------------------------------------------------
# beta <-> M transformations (no minfi needed)
# ---------------------------------------------------------------------------

.champ_beta2m <- function(beta) {
    beta <- as.matrix(beta)
    m <- log2(beta / (1 - beta))
    m[!is.finite(m)] <- NA_real_
    m
}

.champ_m2beta <- function(m) {
    m <- as.matrix(m)
    2^m / (2^m + 1)
}

.champ_logit2 <- function(x) {
    x <- as.matrix(x)
    log2(x / (1 - x))
}

.champ_ilogit2 <- function(x) {
    x <- as.matrix(x)
    2^x / (2^x + 1)
}

# ---------------------------------------------------------------------------
# matrix / input helpers
# ---------------------------------------------------------------------------

.champ_as_beta <- function(beta, name = "beta") {
    if (is.null(beta)) .champ_stop("'", name, "' must not be NULL.")
    if (is.data.frame(beta)) beta <- as.matrix(beta)
    if (!is.matrix(beta)) {
        if (is.vector(beta)) {
            beta <- matrix(beta, ncol = 1L, dimnames = list(names(beta), NULL))
        } else {
            .champ_stop("'", name, "' must be a matrix or data.frame.")
        }
    }
    storage.mode(beta) <- "double"
    if (is.null(rownames(beta))) .champ_stop("'", name, "' must have rownames (probe IDs).")
    if (is.null(colnames(beta))) .champ_stop("'", name, "' must have colnames (sample IDs).")
    beta
}

.champ_cores <- function(cores = 1L) {
    if (!is.numeric(cores) || length(cores) != 1L || is.na(cores)) cores <- 1L
    cores <- as.integer(cores)
    if (cores < 1L) cores <- 1L
    avail <- .champ_try(parallel::detectCores(logical = FALSE), default = 1L)
    if (is.null(avail) || is.na(avail) || avail < 1L) avail <- 1L
    if (cores > avail) {
        .champ_msg("Requested ", cores, " cores but only ", avail, " are available; using ", avail, ".")
        cores <- avail
    }
    cores
}

#' Coerce a phenotype vector to character (factors keep their labels).
.champ_pheno <- function(pheno, n = NULL, name = "pheno") {
    if (is.null(pheno)) .champ_stop("'", name, "' must not be NULL.")
    if (is.factor(pheno)) pheno <- as.character(pheno)
    if (is.data.frame(pheno) || is.matrix(pheno)) {
        if (ncol(pheno) == 1L) {
            pheno <- pheno[[1L]]
        } else {
            .champ_stop("'", name, "' must be a vector, not a ", nrow(pheno), " x ", ncol(pheno), " table.")
        }
    }
    if (!is.null(n) && length(pheno) != n) {
        .champ_stop("'", name, "' has length ", length(pheno), " but the data has ", n, " samples.")
    }
    pheno
}

#' Grouping factor with a deterministic level order.
.champ_group <- function(pheno, levels = NULL) {
    pheno <- .champ_pheno(pheno)
    factor(pheno, levels = if (is.null(levels)) sort(unique(pheno)) else levels)
}

.champ_is_numeric <- function(x) is.numeric(x) && !is.factor(x)

#' Row-wise variance (matrixStats is optional).
.champ_row_vars <- function(x) {
    if (.champ_has("matrixStats")) return(matrixStats::rowVars(x))
    apply(x, 1L, stats::var)
}

#' Column medians / MADs that work across matrixStats versions.
.champ_col_medians <- function(x) {
    if (.champ_has("matrixStats")) return(matrixStats::colMedians(x))
    apply(x, 2L, stats::median)
}

.champ_col_mads <- function(x) {
    if (.champ_has("matrixStats")) {
        f <- .champ_try(getExportedValue("matrixStats", "colMads"), default = NULL)
        if (!is.null(f)) return(f(x))
        f <- .champ_try(getExportedValue("matrixStats", "colMad"), default = NULL)
        if (!is.null(f)) return(f(x))
    }
    apply(x, 2L, stats::mad)
}

#' Safe access to a list element by name, returning NULL instead of erroring.
.champ_elt <- function(x, name) if (is.list(x) && name %in% names(x)) x[[name]] else NULL

# ---------------------------------------------------------------------------
# dependency manifest / environment report
# ---------------------------------------------------------------------------

#' The canonical dependency manifest (kept in sync with DESCRIPTION).
.champ_dependency_manifest <- function() {
    imports <- c(
        "ChAMPdata", "limma", "matrixStats", "illuminaio", "DNAcopy",
        "preprocessCore", "bumphunter", "sva", "globaltest", "GenomicRanges",
        "IRanges", "S4Vectors", "GenomeInfoDb", "quadprog")
    suggests <- c(
        "minfi", "DMRcate", "impute", "RPMM", "kpmt", "goseq",
        "geneLenDataBase", "missMethyl", "dendextend", "shiny", "shinythemes",
        "plotly", "DT", "prettydoc", "rmarkdown", "knitr", "testthat",
        "BiocManager", "org.Hs.eg.db")
    bioc <- c("ChAMPdata", "limma", "DNAcopy", "preprocessCore", "bumphunter",
              "sva", "globaltest", "GenomicRanges", "IRanges", "S4Vectors",
              "GenomeInfoDb", "illuminaio", "minfi", "DMRcate", "impute",
              "goseq", "geneLenDataBase", "missMethyl", "org.Hs.eg.db")
    archived <- c("kpmt")
    fmt <- function(p, role) {
        src <- ifelse(p %in% archived, "CRAN (archived)",
                      ifelse(p %in% bioc, "Bioconductor", "CRAN"))
        data.frame(package = p, role = role, source = src,
                   stringsAsFactors = FALSE)
    }
    rbind(fmt(imports, "Imports"), fmt(suggests, "Suggests"))
}

#' Data frame describing the ChAMP dependency status.
.champ_dependency_table <- function() {
    man <- .champ_dependency_manifest()
    rows <- lapply(seq_len(nrow(man)), function(i) {
        p  <- man$package[i]
        ok <- .champ_has(p)
        data.frame(
            package   = p,
            role      = man$role[i],
            installed = ok,
            version   = if (ok) as.character(utils::packageVersion(p)) else NA_character_,
            source    = man$source[i],
            stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}