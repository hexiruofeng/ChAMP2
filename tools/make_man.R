#!/usr/bin/env Rscript
# =============================================================================
#  make_man.R -- generate man/*.Rd for every exported ChAMP function.
#
#  The package documents its functions with roxygen-style comment blocks
#  (@param/@return), but the Rd files are materialised by this script instead
#  of by roxygen2 so that no extra build dependency is required.  Argument
#  descriptions that have no roxygen entry fall back to the dictionary below.
#
#  Usage (from the package root):
#      Rscript tools/make_man.R
# =============================================================================

options(warn = 1)
root <- normalizePath(".")
r_files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)

## ---------------------------------------------------------------- roxygen ---
parse_blocks <- function(files) {
    out <- list()
    for (f in files) {
        lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
        cur <- character(0)
        for (i in seq_along(lines)) {
            l <- lines[[i]]
            if (grepl("^#'", l)) {
                cur <- c(cur, sub("^#'[ ]?", "", l))
                next
            }
            if (grepl("^[A-Za-z._][A-Za-z0-9._]*[ ]*<-[ ]*function", l) &&
                length(cur)) {
                nm <- sub("[ ]*<-.*$", "", trimws(l))
                out[[nm]] <- cur
                cur <- character(0)
                next
            }
            if (!grepl("^[[:space:]]*$", l)) cur <- character(0)
        }
    }
    out
}

tidy <- function(x) {
    x <- trimws(x)
    x <- x[nzchar(x)]
    paste(x, collapse = " ")
}

parse_block <- function(lines) {
    blank <- which(!nzchar(trimws(lines)))
    end_title <- if (length(blank)) blank[1L] - 1L else length(lines)
    title <- tidy(lines[seq_len(max(0L, end_title))])
    body <- if (length(blank)) lines[(blank[1L] + 1L):length(lines)] else character(0)
    if (!length(body)) body <- character(0)
    tags <- grep("^@", body)
    desc_txt <- if (length(tags)) {
        tidy(body[seq_len(max(0L, tags[1L] - 1L))])
    } else {
        tidy(body)
    }
    res <- list(title = title, description = desc_txt, params = character(0),
                return = "", details = "")
    lines <- body
    tags <- grep("^@", lines)
    if (length(tags)) {
        for (t in tags) {
            body <- lines[t]
            tag <- sub("^@([A-Za-z]+).*$", "\\1", body)
            rest <- trimws(sub("^@[A-Za-z]+[ ]?", "", body))
            cont <- if (t < length(lines)) {
                nxt <- if (length(tags)) tags[tags > t][1] else NA_integer_
                end <- if (is.na(nxt)) length(lines) else nxt - 1L
                if (end > t) lines[(t + 1L):end] else character(0)
            } else character(0)
            cont <- cont[!grepl("^@", cont)]
            txt <- tidy(c(rest, cont))
            if (identical(tag, "param")) {
                parts <- strsplit(txt, "[ ]+")[[1]]
                nm <- parts[1]
                desc <- tidy(paste(parts[-1], collapse = " "))
                res$params[nm] <- desc
            } else if (identical(tag, "return")) {
                res$return <- txt
            } else if (identical(tag, "description")) {
                res$description <- txt
            } else if (identical(tag, "details")) {
                res$details <- txt
            }
        }
    }
    res
}

## ---------------------------------------------------------------- exports ---
exported <- local({
    ns <- readLines(file.path(root, "NAMESPACE"))
    txt <- paste(ns, collapse = "\n")
    m <- regmatches(txt, regexpr("export\\([^)]*\\)", txt))
    if (!length(m)) return(character(0))
    inner <- sub("^export\\(([^)]*)\\)$", "\\1", m)
    x <- trimws(strsplit(inner, ",")[[1]])
    x[nzchar(x)]
})

## --------------------------------------------------------------- dictionary -
arg_dict <- c(
    beta = "A numeric matrix of beta values, with probes in rows and samples in columns.",
    beta.m = "A numeric matrix of beta values (probes x samples).",
    beta.v = "Beta values of a single sample (NAs are not allowed).",
    M = "A numeric matrix of M values (probes x samples).",
    pd = "A data.frame with one row per sample (the sample sheet).",
    pheno = "A phenotype vector, one value per sample. Character, factor or numeric.",
    arraytype = "Array type: '450K', 'EPICv1' or 'EPICv2' ('EPIC' is an alias of 'EPICv2').",
    verbose = "Print progress messages.",
    method = "The method to use; see the description for the accepted values.",
    resultsDir = "Directory where result files are written.",
    cores = "Number of workers to use.",
    intensity = "A matrix of total intensities (meth + unmeth), probes x samples.",
    Meth = "Optional matrix of methylated intensities, filtered alongside the data.",
    UnMeth = "Optional matrix of unmethylated intensities, filtered alongside the data.",
    detP = "A detection p-value matrix with the same dimensions as the data matrix.",
    beadcount = "A bead count matrix with the same dimensions as the data matrix.",
    autoimpute = "Impute remaining missing values with KNN (k = 5).",
    filterDetP = "Filter probes and samples on the detection p-value.",
    ProbeCutoff = "Fraction of samples in which a probe may fail before it is removed (0 = any failure removes it).",
    SampleCutoff = "Fraction of failed probes above which a sample is removed.",
    detPcut = "Detection p-value threshold.",
    filterBeads = "Filter probes with bead counts below three.",
    beadCutoff = "Fraction of samples with bead count < 3 above which a probe is removed.",
    filterNoCG = "Remove non-CpG probes.",
    filterSNPs = "Remove probes overlapping known SNPs.",
    population = "Population used for the SNP mask (e.g. 'EUR', 'AFR'); NULL uses the general mask.",
    filterMultiHit = "Remove probes that map to multiple genomic locations.",
    filterXY = "Remove probes on the X and Y chromosomes.",
    fixOutlier = "Replace exact 0 and 1 beta values by the closest valid value.",
    directory = "Directory containing the IDAT files (and a sample sheet for method 'minfi').",
    methValue = "'B' for beta values or 'M' for M values.",
    force = "Passed to minfi::read.metharray.exp(); read incomplete arrays.",
    k = "Number of nearest neighbours.",
    imputemethod = "Imputation method: 'Combine', 'Delete' or 'KNN'.",
    RGEffect = "Include control probe intensities as covariates in the SVD.",
    batchname = "Character vector of batch variables to correct.",
    variablename = "Name of the phenotype column that must be protected.",
    logitTrans = "Apply ComBat on the logit scale (recommended for beta values).",
    compare.group = "Two group labels to compare; NULL compares every pair.",
    adjPVal = "Adjusted p-value cutoff used when reporting results.",
    adjust.method = "Multiple testing correction passed to limma::topTable().",
    minProbes = "Minimum number of probes per region.",
    adjPvalDmr = "Adjusted p-value cutoff for reported regions.",
    maxGap = "Maximum distance (bp) between probes in the same cluster.",
    cutoff = "bumphunter cutoff; NULL lets bumphunter pick one.",
    pickCutoff = "Let bumphunter choose the cutoff.",
    smooth = "Smooth the test statistics along the genome.",
    smoothFunction = "Smoothing function passed to bumphunter.",
    useWeights = "Use weights in bumphunter.",
    permutations = "Pre-computed permutations for bumphunter.",
    B = "Number of bumphunter bootstrap permutations.",
    nullMethod = "Null distribution method for bumphunter.",
    meanLassoRadius = "Target mean lasso radius (bp) for ProbeLasso.",
    minDmrSep = "Minimum gap (bp) between merged ProbeLasso regions.",
    minDmrSize = "Minimum DMR width (bp) for ProbeLasso.",
    adjPvalProbe = "Probe-level adjusted p-value cutoff for ProbeLasso.",
    Rplot = "Draw plots on the current device.",
    PDFplot = "Write plots to PDF files in resultsDir.",
    rmSNPCH = "Remove SNP and cross-hybridising probes before DMRcate.",
    dist = "Distance (bp) used when removing SNP probes for DMRcate.",
    mafcut = "Minor allele frequency cutoff for SNP probe removal.",
    fdr = "Probe-level FDR used by DMRcate.",
    lambda = "Kernel bandwidth (bp) for DMRcate.",
    C = "Bandwidth scaling for DMRcate (sigma = lambda / C).",
    dmrcate.backend = "'auto', 'package' or 'internal': which DMRcate implementation to use.",
    DMP = "A champ.DMP() result (a data.frame or a named list of data.frames).",
    DMR = "A champ.DMR() result (a data.frame or a named list of data.frames).",
    CpGlist = "Character vector or named list of CpG identifiers.",
    Genelist = "Character vector or named list of gene symbols.",
    adjPval = "Adjusted p-value cutoff.",
    min.cpgs = "Minimum number of probes per DMR.",
    pcutoff = "'fdr' or a numeric probe-level p-value cutoff.",
    betacutoff = "Optional absolute mean beta difference filter.",
    control = "Use control samples to compute the CNA log ratios.",
    controlGroup = "Name of the control group; 'champCtls' uses the bundled blood controls.",
    sampleCNA = "Segment every sample individually.",
    groupFreqPlots = "Draw gain/loss frequency plots per group.",
    freqThreshold = "Frequency threshold for the group CNA plots.",
    plotBMIQ = "Write BMIQ diagnostic PDFs.",
    mset = "A minfi MethylSet (required for SWAN).",
    rgSet = "A minfi RGChannelSet (required for SWAN and FunctionalNormalization).",
    nL = "Number of mixture components (BMIQ requires 3).",
    doH = "Correct hemimethylated type II probes.",
    nfit = "Number of probes used to fit the beta mixture.",
    th1.v = "Initial thresholds for type I probes.",
    th2.v = "Initial thresholds for type II probes (estimated when NULL).",
    niter = "Number of EM iterations.",
    tol = "Convergence tolerance.",
    plots = "Write diagnostic plots.",
    plotDir = "Directory for diagnostic plots.",
    sampleID = "Sample label used in plot file names.",
    design.v = "Probe design vector: 1 = type I, 2 = type II.",
    pnbeta.v = "BMIQ normalised beta values.",
    Y = "A numeric matrix (n x J), typically one column of beta values.",
    w = "An n x K matrix of initial posterior class memberships.",
    maxiter = "Maximum number of EM iterations.",
    weights = "Case weights.",
    pop = "Named numeric vector of population values.",
    obs = "Character vector of population member names that were observed.",
    med = "Pre-computed minimal median (alternative to obs).",
    size = "Sample size; required when med is supplied.",
    tail = "One of 'lower', 'upper' or 'two-sided'.",
    data.m = "A numeric matrix (probes x samples).",
    betaFile = "Tab separated beta matrix file.",
    sampleSheet = "Tab separated sample sheet.",
    sep = "Field separator.",
    dat = "A data matrix (features x samples).",
    batch = "A batch vector, one value per sample.",
    mod = "A model matrix of covariates (a 'Batch' column is removed).",
    numCovs = "Unused; kept for backwards compatibility.",
    par.prior = "Use the parametric empirical Bayes prior.",
    prior.plots = "Draw prior diagnostic plots.",
    mean.only = "Only correct the mean (sva >= 3.30).",
    optional = "Also install the Suggests packages.",
    binary = "Use binary packages when available.",
    upgrade = "Passed to install.packages(): 'never', 'default' or 'always'.",
    force_reinstall = "Reinstall packages that are already present.",
    which = "'Imports', 'Suggests' or 'all'.",
    quiet = "Do not print the report.",
    runApp = "Launch the Shiny application; FALSE returns the app object instead.",
    cutgroupnumber = "Number of quantile groups used when pheno is numeric.",
    CpG = "Character vector of CpG identifiers to summarise.",
    Block = "A champ.Block() result.")

## --------------------------------------------------------------- helpers ----
rd_escape <- function(x) {
    x <- gsub("\\\\", "\\\\\\\\", x)
    x <- gsub("%", "\\\\%", x)
    x <- gsub("\\{", "\\\\{", x)
    x <- gsub("\\}", "\\\\}", x)
    x
}

fmt_usage <- function(name, fn) {
    txt <- paste(utils::capture.output(base::args(fn)), collapse = " ")
    txt <- sub("^function[ ]*", "", txt)
    txt <- sub("[ ]*NULL[ ]*$", "", txt)
    u <- paste0(name, txt)
    u <- gsub("[ ]+", " ", u)
    if (nchar(u) > 78) {
        u <- paste(strwrap(u, width = 76, exdent = 4L), collapse = "\n")
    }
    rd_escape(u)
}

make_rd <- function(name, fn, block) {
    params <- block$params
    f <- formals(fn)
    args <- names(f)
    lines <- c(
        paste0("\\name{", name, "}"),
        paste0("\\alias{", name, "}"),
        paste0("\\title{", rd_escape(if (nzchar(block$title)) block$title else name), "}"),
        "\\usage{",
        fmt_usage(name, fn),
        "}",
        "\\arguments{")
    for (a in args) {
        desc <- if (!is.na(params[a]) && nzchar(params[a])) params[a] else
            if (a %in% names(arg_dict)) arg_dict[[a]] else
                "See the package vignette for details."
        if (identical(a, "...")) desc <- "Further arguments passed to the underlying methods."
        lines <- c(lines, paste0("  \\item{", rd_escape(a), "}{", rd_escape(desc), "}"))
    }
    lines <- c(lines, "}")
    desc <- block$description
    if (!nzchar(desc)) desc <- block$title
    if (!nzchar(desc)) desc <- paste0("See the ChAMP vignette for a description of ", name, ".")
    lines <- c(lines, paste0("\\description{", rd_escape(desc), "}"))
    if (nzchar(block$details)) {
        lines <- c(lines, paste0("\\details{", rd_escape(block$details), "}"))
    }
    ret <- if (nzchar(block$return)) block$return else
        "See the description above; the returned object is documented in the vignette."
    lines <- c(lines,
               paste0("\\value{", rd_escape(ret), "}"),
               "\\seealso{\\code{\\link{champ.process}}, \\code{\\link{champ_check_env}}}",
               "\\keyword{methods}")
    lines
}

## ------------------------------------------------------------------- main ---
blocks <- lapply(parse_blocks(r_files), parse_block)
env <- new.env(parent = globalenv())
for (f in r_files) sys.source(f, envir = env)

man_dir <- file.path(root, "man")
dir.create(man_dir, showWarnings = FALSE)

written <- character(0)
for (nm in exported) {
    if (!exists(nm, envir = env, inherits = FALSE)) {
        message("skip (not found): ", nm)
        next
    }
    fn <- get(nm, envir = env)
    if (!is.function(fn)) next
    blk <- if (!is.null(blocks[[nm]])) blocks[[nm]] else
        list(title = nm, description = "", params = character(0), return = "", details = "")
    writeLines(make_rd(nm, fn, blk), file.path(man_dir, paste0(nm, ".Rd")),
               useBytes = TRUE)
    written <- c(written, nm)
}

## package-level page
pkg <- c(
    "\\name{ChAMP-package}",
    "\\alias{ChAMP-package}",
    "\\alias{ChAMP}",
    "\\docType{package}",
    "\\title{ChAMP2: Chip Analysis Methylation Pipeline}",
    "\\description{",
    "  Quality control, normalisation, differential methylation, copy number and",
    "  gene set enrichment analysis for Illumina HumanMethylation 450K, EPIC and",
    "  EPICv2 arrays.  This fork (ChAMP2) installs and runs on R 4.2.3 to 4.6.2.",
    "}",
    "\\details{",
    "  The main entry point is \\code{\\link{champ.process}}; individual steps can be",
    "  run with \\code{\\link{champ.load}}, \\code{\\link{champ.filter}},",
    "  \\code{\\link{champ.QC}}, \\code{\\link{champ.norm}}, \\code{\\link{champ.SVD}},",
    "  \\code{\\link{champ.runCombat}}, \\code{\\link{champ.DMP}}, \\code{\\link{champ.DMR}},",
    "  \\code{\\link{champ.Block}}, \\code{\\link{champ.GSEA}}, \\code{\\link{champ.ebGSEA}},",
    "  \\code{\\link{champ.CNA}} and \\code{\\link{champ.refbase}}.",
    "  Run \\code{\\link{champ_check_env}} to see which optional features are available",
    "  and \\code{\\link{champ_install_deps}} to install the missing packages.",
    "}",
    "\\keyword{package}")
writeLines(pkg, file.path(man_dir, "ChAMP-package.Rd"), useBytes = TRUE)

cat("Generated", length(written), "Rd files in", man_dir, "\n")