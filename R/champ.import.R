# =============================================================================
#  champ.import -- read Illumina IDAT files with the "ChAMP" (annotation based)
#  method.  No minfi involvement, which makes this path usable on every
#  supported R / Bioconductor combination.
# =============================================================================

#' Import Illumina IDAT files using the ChAMP annotation
#'
#' Reads a sample sheet and the matching Grn/Red IDAT files, and returns beta
#' values, M values, intensities, detection p-values and bead counts.
#'
#' @param directory directory containing the sample sheet (csv) and IDAT files.
#' @param offset offset added to the denominator when computing beta values.
#' @param arraytype "450K", "EPICv1" or "EPICv2" ("EPIC" is an alias of
#'   "EPICv2").
#' @param verbose print progress.
#'
#' @return list with beta, M, pd, intensity, detP, beadcount, Meth, UnMeth,
#'   arraytype.
#' @export
champ.import <- function(directory = getwd(), offset = 100, arraytype = "450K",
                         verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<< ChAMP.IMPORT START >>>>>]")
    .champ_msg("-----------------------------")
    at <- .champ_arraytype(arraytype)
    .champ_require("illuminaio", "champ.import()")

    .champ_msg("\n[ Section 1: Read PD Files Start ]")
    if (!dir.exists(directory)) {
        .champ_stop("  Your 'directory' for loading does not exist: ", directory)
    }
    csvfile <- list.files(directory, recursive = TRUE, pattern = "[.]csv$",
                          full.names = TRUE, ignore.case = TRUE)
    if (length(csvfile) == 0L) {
        .champ_stop("  champ.import can not find any csv file in ", directory, ".")
    } else if (length(csvfile) > 1L) {
        .champ_stop("  champ.import finds more than one csv file in ", directory, ".")
    }
    .champ_msg("  CSV Directory: ", csvfile)

    lines <- readLines(csvfile, warn = FALSE)
    skipline <- which(substr(trimws(lines), 1L, 6L) == "[Data]")
    pd <- suppressWarnings(utils::read.csv(csvfile, skip = skipline,
                                           header = TRUE,
                                           stringsAsFactors = FALSE))
    if (!is.data.frame(pd) || nrow(pd) == 0L) {
        .champ_stop("  The sample sheet ", csvfile, " contains no data.")
    }

    if ("Sentrix_Position" %in% colnames(pd)) {
        colnames(pd)[colnames(pd) == "Sentrix_Position"] <- "Array"
        .champ_msg("  Replace Sentrix_Position into Array")
    }
    if ("Sentrix_ID" %in% colnames(pd)) {
        colnames(pd)[colnames(pd) == "Sentrix_ID"] <- "Slide"
        pd$Slide <- as.character(pd$Slide)
        .champ_msg("  Replace Sentrix_ID into Slide")
    }
    for (nm in c("Pool_ID", "Sample_Plate", "Sample_Well", "Sample_Name", "Sample_Group")) {
        if (nm %in% colnames(pd)) pd[[nm]] <- as.character(pd[[nm]])
    }
    if (!all(c("Slide", "Array") %in% colnames(pd))) {
        .champ_stop("  Your sample sheet must contain Slide (Sentrix_ID) and Array (Sentrix_Position).")
    }
    if (!"Sample_Name" %in% colnames(pd)) {
        pd$Sample_Name <- paste(pd$Slide, pd$Array, sep = "_")
    }

    .champ_msg("[ Section 1: Read PD file Done ]")
    .champ_msg("\n[ Section 2: Read IDAT files Start ]")

    all_files <- list.files(directory, recursive = TRUE, full.names = TRUE,
                            pattern = "[.]idat$", ignore.case = TRUE)
    grn_idx <- grep("_Grn[.]idat$", all_files, ignore.case = TRUE)
    red_idx <- grep("_Red[.]idat$", all_files, ignore.case = TRUE)
    grn_files <- all_files[grn_idx]
    red_files <- all_files[red_idx]
    names(grn_files) <- sub("_Grn[.]idat$", "", basename(grn_files), ignore.case = TRUE)
    names(red_files) <- sub("_Red[.]idat$", "", basename(red_files), ignore.case = TRUE)

    key <- paste(pd$Slide, pd$Array, sep = "_")
    GrnPath <- grn_files[key]
    RedPath <- red_files[key]
    if (anyNA(GrnPath) || anyNA(RedPath)) {
        missing <- key[is.na(GrnPath) | is.na(RedPath)]
        .champ_stop("  No matching IDAT files were found for: ",
                    paste(utils::head(missing, 5L), collapse = ", "),
                    "\n  Found ", length(grn_files), " Grn and ",
                    length(red_files), " Red IDAT files.")
    }

    read_one <- function(path, tag) {
        .champ_msg("  Loading ", tag, ": ", basename(path))
        illuminaio::readIDAT(path)
    }
    G.idats <- lapply(GrnPath, read_one, tag = "Grn")
    R.idats <- lapply(RedPath, read_one, tag = "Red")
    names(G.idats) <- pd$Sample_Name
    names(R.idats) <- pd$Sample_Name

    commonAddresses <- as.character(Reduce(intersect,
        lapply(c(G.idats, R.idats), function(x) rownames(x$Quants))))
    if (!length(commonAddresses)) {
        .champ_stop("  The IDAT files have no addresses in common.")
    }

    GreenMean <- do.call(cbind, lapply(G.idats, function(x)
        x$Quants[commonAddresses, "Mean"]))
    RedMean <- do.call(cbind, lapply(R.idats, function(x)
        x$Quants[commonAddresses, "Mean"]))
    .champ_msg("    Your Red Green Channel contains ", nrow(GreenMean), " probes.")
    .champ_msg("[ Section 2: Read IDAT Files Done ]")

    .champ_msg("\n[ Section 3: Use Annotation Start ]")
    .champ_msg("\n  Reading ", at, " Annotation >>")
    anno_name <- switch(at,
        "450K"   = "Anno450K",
        "EPICv1" = "AnnoEPICv1",
        "EPICv2" = "AnnoEPICv2",
        "27K"    = "Anno450K",
        .champ_stop("champ.import() supports 450K, EPICv1 and EPICv2."))
    Anno <- .champ_data_object(anno_name)
    if (!is.list(Anno) || is.null(Anno$Annotation)) {
        .champ_stop("  ", anno_name, " does not contain an $Annotation table.")
    }
    annotation <- as.data.frame(Anno$Annotation, stringsAsFactors = FALSE)
    cpg_col <- intersect(c("CpG", "Name", "IlmnID", "probeID"), colnames(annotation))
    if (length(cpg_col)) {
        annotation$CpG <- as.character(annotation[[cpg_col[1L]]])
    } else if (!is.null(rownames(Anno$Annotation))) {
        annotation$CpG <- rownames(Anno$Annotation)
    } else {
        .champ_stop("  ", anno_name, " has no CpG identifier column.")
    }
    for (nm in c("M.index", "U.index", "Channel")) {
        if (!nm %in% colnames(annotation)) {
            .champ_stop("  ", anno_name, " is missing the '", nm, "' column.")
        }
    }

    control_probe <- character(0)
    if (!is.null(Anno$ControlProbe)) {
        cp <- as.data.frame(Anno$ControlProbe, stringsAsFactors = FALSE)
        control_probe <- rownames(cp)[cp[, 1L] == "NEGATIVE"]
        control_probe <- control_probe[control_probe %in% rownames(RedMean)]
        .champ_msg("    ", length(control_probe), " negative control probes found.")
    }

    R.Load <- RedMean
    G.Load <- GreenMean
    rownames(G.Load) <- paste0("G-", rownames(G.Load))
    rownames(R.Load) <- paste0("R-", rownames(R.Load))
    IDAT <- rbind(G.Load, R.Load)

    if (!length(control_probe)) {
        ## fall back to a robust background estimate
        rMu <- rep(stats::median(RedMean), ncol(RedMean))
        rSd <- rep(stats::mad(as.vector(RedMean)), ncol(RedMean))
        gMu <- rep(stats::median(GreenMean), ncol(GreenMean))
        gSd <- rep(stats::mad(as.vector(GreenMean)), ncol(GreenMean))
    } else {
        rMu <- .champ_col_medians(R.Load[control_probe, , drop = FALSE])
        rSd <- .champ_col_mads(R.Load[control_probe, , drop = FALSE])
        gMu <- .champ_col_medians(G.Load[control_probe, , drop = FALSE])
        gSd <- .champ_col_mads(G.Load[control_probe, , drop = FALSE])
    }
    rSd[!is.finite(rSd) | rSd <= 0] <- 1
    gSd[!is.finite(gSd) | gSd <= 0] <- 1

    .champ_msg("\n  Generating Meth and UnMeth Matrix")
    m_ok <- annotation$M.index %in% rownames(IDAT)
    u_ok <- annotation$U.index %in% rownames(IDAT)
    if (!identical(m_ok, u_ok)) {
        .champ_stop("  Meth and UnMeth matrices are not paired correctly.")
    }
    if (!any(m_ok)) .champ_stop("  None of the annotation addresses were found in the IDAT files.")
    M <- IDAT[annotation$M.index[m_ok], , drop = FALSE]
    U <- IDAT[annotation$U.index[m_ok], , drop = FALSE]
    cpg_index <- annotation$CpG[m_ok]
    keep <- !duplicated(cpg_index)
    M <- M[keep, , drop = FALSE]
    U <- U[keep, , drop = FALSE]
    rownames(M) <- rownames(U) <- cpg_index[keep]

    .champ_msg("\n  Generating beta Matrix")
    BetaValue <- M / (M + U + offset)
    .champ_msg("  Generating M Matrix")
    MValue <- log2(M / U)
    .champ_msg("  Generating intensity Matrix")
    intensity <- M + U

    .champ_msg("  Calculating Detect P value")
    detP <- matrix(NA_real_, nrow = nrow(intensity), ncol = ncol(intensity),
                   dimnames = dimnames(intensity))
    channel <- annotation$Channel[m_ok][keep]
    type_II   <- rownames(detP)[channel == "g+r"]
    type_I_g  <- rownames(detP)[channel %in% c("g", "Grn", "G")]
    type_I_r  <- rownames(detP)[channel %in% c("r", "Red", "R")]
    for (i in seq_len(ncol(detP))) {
        if (length(type_II)) {
            detP[type_II, i] <- 1 - stats::pnorm(intensity[type_II, i],
                                                 mean = rMu[i] + gMu[i],
                                                 sd = rSd[i] + gSd[i])
        }
        if (length(type_I_r)) {
            detP[type_I_r, i] <- 1 - stats::pnorm(intensity[type_I_r, i],
                                                  mean = rMu[i] * 2, sd = rSd[i] * 2)
        }
        if (length(type_I_g)) {
            detP[type_I_g, i] <- 1 - stats::pnorm(intensity[type_I_g, i],
                                                  mean = gMu[i] * 2, sd = gSd[i] * 2)
        }
    }
    detP[detP < 0] <- 0
    if (anyNA(detP)) .champ_msg("    !!! There are NA values in your detP matrix.")

    .champ_msg("  Counting Beads")
    NBeads <- do.call(cbind, lapply(R.idats, function(x)
        x$Quants[commonAddresses, "NBeads"]))
    m_addr <- sub("^[GR]-", "", annotation$M.index[m_ok][keep])
    u_addr <- sub("^[GR]-", "", annotation$U.index[m_ok][keep])
    Mbead <- NBeads[m_addr, , drop = FALSE]
    Ubead <- NBeads[u_addr, , drop = FALSE]
    Ubead[Ubead < 3 | Mbead < 3] <- NA_real_
    rownames(Ubead) <- rownames(intensity)

    .champ_msg("[ Section 3: Use Annotation Done ]")
    .champ_msg("\n[<<<<< ChAMP.IMPORT END >>>>>>]")
    .champ_msg("[===========================]")
    .champ_msg("[You may want to process champ.filter() next.]\n")

    list(beta = BetaValue, M = MValue, pd = pd, intensity = intensity,
         detP = detP, beadcount = Ubead, Meth = M, UnMeth = U,
         arraytype = at)
}