# =============================================================================
#  champ.GSEA -- gene set enrichment analysis
# =============================================================================

#' Collapse a vector of gene symbols into one string
#' @keywords internal
.champ_paste_vector <- function(v) {
    v <- v[!is.na(v) & nzchar(v)]
    if (!length(v)) return("")
    paste(unique(v), collapse = " ")
}

#' Gene set enrichment analysis
#'
#' Four methods are supported:
#' \describe{
#'   \item{fisher}{Fisher's exact test on MSigDB-like gene sets, weighted by
#'     the number of CpGs per gene (fast, no extra packages).}
#'   \item{goseq}{goseq, which models the CpG-count bias explicitly.}
#'   \item{gometh}{missMethyl::gometh for GO enrichment of CpG lists.}
#'   \item{ebayes}{champ.ebGSEA(), a global-test based enrichment analysis.}
#' }
#'
#' @param beta beta value matrix (used for annotation and, for "ebayes", the
#'   test itself).
#' @param DMP a champ.DMP() result (single data.frame) or NULL.
#' @param DMR a champ.DMR() result or NULL.
#' @param CpGlist character vector or named list of CpG identifiers.
#' @param Genelist character vector or named list of gene symbols.
#' @param pheno phenotype vector (required for method "ebayes").
#' @param method "fisher", "goseq", "gometh" or "ebayes".
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param Rplot draw diagnostic plots.
#' @param adjPval adjusted p-value cutoff.
#' @param cores number of workers for the "ebayes" method.
#' @param verbose print progress.
#'
#' @return named list of enrichment tables.
#' @export
champ.GSEA <- function(beta = NULL, DMP = NULL, DMR = NULL, CpGlist = NULL,
                       Genelist = NULL, pheno = NULL, method = "fisher",
                       arraytype = "450K", Rplot = FALSE, adjPval = 0.05,
                       cores = 1L, verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<< ChAMP.GSEA START >>>>>]")
    .champ_msg("-----------------------------")
    beta <- .champ_as_beta(beta)
    at <- .champ_arraytype(arraytype)
    method <- match.arg(method, c("fisher", "goseq", "gometh", "ebayes"))
    pf <- .champ_probe_features(at)
    RSanno <- pf[rownames(beta), c("CHR", "MAPINFO", "gene"), drop = FALSE]
    colnames(RSanno) <- c("chr", "pos", "UCSC_RefGene_Name")
    RSanno$UCSC_RefGene_Name[is.na(RSanno$UCSC_RefGene_Name)] <- ""

    .champ_require("stats", "champ.GSEA()")
    if (identical(method, "ebayes")) {
        if (is.null(pheno)) .champ_stop("'pheno' is required for method = 'ebayes'.")
        out <- champ.ebGSEA(beta = beta, pheno = pheno, adjPval = adjPval,
                            arraytype = at, cores = cores, verbose = verbose)
        .champ_msg("[<<<<< ChAMP.GSEA END >>>>>>]")
        return(out)
    }

    if (identical(method, "gometh")) {
        .champ_require("missMethyl", "method = 'gometh'")
        genes <- NULL
        loi.lv <- list()
        if (!is.null(DMP)) {
            DMP <- if (is.list(DMP) && !is.data.frame(DMP)) DMP[[1L]] else DMP
            loi.lv[["DMP"]] <- rownames(DMP)
        }
        if (!is.null(DMR)) {
            dmr <- if (is.data.frame(DMR)) DMR else DMR[[1L]]
            loi.lv[["DMR"]] <- .champ_cpgs_in_regions(dmr, RSanno)
        }
        if (!is.null(CpGlist)) loi.lv[["CpG"]] <- as.character(CpGlist)
        if (!is.null(Genelist)) {
            .champ_msg("gometh() only accepts CpGs; 'Genelist' is ignored. Use method = 'fisher'.")
        }
        if (!length(loi.lv)) .champ_stop("Nothing to test: supply DMP, DMR or CpGlist.")
        out <- list()
        for (nm in names(loi.lv)) {
            .champ_msg("  gometh on ", nm, " (", length(loi.lv[[nm]]), " CpGs)")
            tab <- .champ_try(missMethyl::gometh(
                sig.cpg = loi.lv[[nm]], all.cpg = rownames(beta),
                collection = c("GO"), array.type = at, plot.bias = Rplot,
                prior.prob = TRUE), default = NULL)
            if (is.null(tab)) {
                .champ_warn("gometh failed for ", nm, "; skipping.")
                next
            }
            tab <- tab[tab$FDR <= adjPval, , drop = FALSE]
            out[[nm]] <- tab[order(tab$FDR), , drop = FALSE]
        }
        .champ_msg("[<<<<< ChAMP.GSEA END >>>>>>]")
        return(out)
    }

    PathwayList <- .champ_pathway_list()
    ueid.v <- unique(unlist(strsplit(RSanno$UCSC_RefGene_Name, ";", fixed = TRUE),
                            use.names = FALSE))
    ueid.v <- ueid.v[nzchar(ueid.v) & !is.na(ueid.v)]
    bias.Data <- table(unlist(strsplit(RSanno$UCSC_RefGene_Name, ";", fixed = TRUE),
                              use.names = FALSE))
    bias.Data <- as.numeric(bias.Data[ueid.v])

    loi.lv <- list()
    if (!is.null(DMP)) {
        DMP <- if (is.list(DMP) && !is.data.frame(DMP)) DMP[[1L]] else DMP
        cpg.idx <- rownames(DMP)[DMP$adj.P.Val <= adjPval]
        loi.lv[["DMP"]] <- .champ_genes_of_cpgs(cpg.idx, RSanno)
    }
    if (!is.null(DMR)) {
        dmr <- if (is.data.frame(DMR)) DMR else DMR[[1L]]
        loi.lv[["DMR"]] <- .champ_genes_of_cpgs(.champ_cpgs_in_regions(dmr, RSanno), RSanno)
    }
    if (!is.null(CpGlist)) {
        if (is.character(CpGlist)) {
            loi.lv[["CpG"]] <- .champ_genes_of_cpgs(CpGlist, RSanno)
        } else if (is.list(CpGlist)) {
            if (is.null(names(CpGlist)) || any(is.na(names(CpGlist)))) {
                .champ_stop("Please name every element of 'CpGlist'.")
            }
            bad <- intersect(names(CpGlist), c("DMP", "DMR", "CpG", "Gene"))
            if (length(bad)) .champ_stop("'CpGlist' must not contain names: ", paste(bad, collapse = ", "))
            loi.lv <- c(loi.lv, lapply(CpGlist, function(x) .champ_genes_of_cpgs(as.character(x), RSanno)))
        } else {
            .champ_stop("'CpGlist' must be a character vector or a named list.")
        }
    }
    if (!is.null(Genelist)) {
        if (is.character(Genelist)) {
            loi.lv[["Gene"]] <- Genelist
        } else if (is.list(Genelist)) {
            if (is.null(names(Genelist)) || any(is.na(names(Genelist)))) {
                .champ_stop("Please name every element of 'Genelist'.")
            }
            bad <- intersect(names(Genelist), c("DMP", "DMR", "CpG", "Gene"))
            if (length(bad)) .champ_stop("'Genelist' must not contain names: ", paste(bad, collapse = ", "))
            loi.lv <- c(loi.lv, Genelist)
        } else {
            .champ_stop("'Genelist' must be a character vector or a named list.")
        }
    }
    if (!length(loi.lv)) .champ_stop("Nothing to test: supply DMP, DMR, CpGlist or Genelist.")
    .champ_msg("<< Gene lists ready: ", paste(names(loi.lv), collapse = ", "), " >>")

    nrepGEID.m <- matrix(NA_real_, nrow = length(PathwayList), ncol = 3L,
                         dimnames = list(names(PathwayList),
                                         c("nList", "nRep", "fRep")))
    nrepGEID.m[, "nList"] <- vapply(PathwayList, length, integer(1L))
    nrepGEID.m[, "nRep"] <- vapply(PathwayList, function(x) length(intersect(x, ueid.v)), integer(1L))
    nrepGEID.m[, "fRep"] <- nrepGEID.m[, "nRep"] / nrepGEID.m[, "nList"]

    listsummary.lm <- list()
    for (i in seq_along(loi.lv)) {
        nm <- names(loi.lv)[i]
        .champ_msg("<< Do GSEA on gene list ", nm, " >>")
        selEID.v <- unique(loi.lv[[nm]])
        selEID.v <- selEID.v[!is.na(selEID.v)]
        if (identical(method, "goseq")) {
            .champ_require("goseq", "method = 'goseq'")
            .champ_require("geneLenDataBase", "method = 'goseq'")
            DE.genes <- as.integer(ueid.v %in% selEID.v)
            names(DE.genes) <- ueid.v
            pwf <- .champ_try(goseq::nullp(DE.genes, "hg19", "geneSymbol",
                                            bias.data = bias.Data, plot.fit = Rplot),
                              default = NULL)
            if (is.null(pwf)) {
                .champ_warn("goseq::nullp() failed for ", nm, "; skipping.")
                next
            }
            pvals <- .champ_try(goseq::goseq(pwf, "hg19", "geneSymbol"), default = NULL)
            if (is.null(pvals)) {
                .champ_warn("goseq() failed for ", nm, "; skipping.")
                next
            }
            pvals$over_represented_adjPvalue <- stats::p.adjust(pvals$over_represented_pvalue,
                                                                method = "BH")
            listsummary.lm[[nm]] <- pvals[pvals$over_represented_adjPvalue <= adjPval, ,
                                          drop = FALSE]
            next
        }

        .champ_msg("  Fisher's exact test (MSigDB gene sets)")
        InMSigDBlist_InList <- vapply(PathwayList, function(x) length(intersect(x, selEID.v)), integer(1L))
        NotInMSigDBlist_InList <- length(selEID.v) - InMSigDBlist_InList
        InMSigDBlist_NotInList <- vapply(PathwayList, function(x)
            length(intersect(x, ueid.v)) - length(intersect(x, selEID.v)), integer(1L))
        NotInMSigDBlist_NotInList <- length(ueid.v) - InMSigDBlist_InList -
            NotInMSigDBlist_InList - InMSigDBlist_NotInList
        ovlapG.lv <- lapply(PathwayList, function(x) intersect(x, selEID.v))
        fisher.m <- cbind(InMSigDBlist_InList, NotInMSigDBlist_InList,
                          InMSigDBlist_NotInList, NotInMSigDBlist_NotInList)
        pv.m <- t(apply(fisher.m, 1L, function(x) {
            ft <- .champ_try(stats::fisher.test(matrix(x, 2L, 2L), alternative = "greater"),
                             default = NULL)
            if (is.null(ft)) c(NA_real_, NA_real_) else c(ft$p.value, unname(ft$estimate))
        }))
        selL.idx <- which(nrepGEID.m[, "fRep"] > 0.6)
        if (!length(selL.idx)) {
            .champ_stop("No gene set covers more than 60% of the array genes. ",
                        "This happens when 'beta' only contains a small subset of the array ",
                        "(use the full matrix) or when the annotation has no gene symbols.")
        }
        out <- data.frame(
            Gene_List = rownames(nrepGEID.m)[selL.idx],
            nrepGEID.m[selL.idx, , drop = FALSE],
            nOVLAP = fisher.m[selL.idx, 1L],
            OR = pv.m[selL.idx, 2L],
            P.value = pv.m[selL.idx, 1L],
            adjPval = stats::p.adjust(pv.m[selL.idx, 1L], method = "BH"),
            Genes = vapply(ovlapG.lv[selL.idx], .champ_paste_vector, character(1L)),
            check.names = FALSE, stringsAsFactors = FALSE)
        out <- out[order(out$adjPval), , drop = FALSE]
        listsummary.lm[[nm]] <- out[out$adjPval <= adjPval, , drop = FALSE]
    }
    .champ_msg("[<<<<< ChAMP.GSEA END >>>>>>]")
    .champ_msg("[===========================]")
    listsummary.lm
}

#' Genes of a set of CpGs
#' @keywords internal
.champ_genes_of_cpgs <- function(cpgs, RSanno) {
    cpgs <- intersect(as.character(cpgs), rownames(RSanno))
    if (!length(cpgs)) return(character(0))
    g <- unlist(strsplit(RSanno[cpgs, "UCSC_RefGene_Name"], ";", fixed = TRUE),
                use.names = FALSE)
    unique(g[nzchar(g) & !is.na(g)])
}

#' CpGs overlapping a DMR table
#' @keywords internal
.champ_cpgs_in_regions <- function(dmr, RSanno) {
    if (is.null(dmr) || !nrow(dmr)) return(character(0))
    if (!all(c("seqnames", "start", "end") %in% colnames(dmr))) {
        .champ_stop("The DMR table must contain seqnames, start and end columns.")
    }
    out <- character(0)
    for (i in seq_len(nrow(dmr))) {
        sel <- RSanno$chr == sub("^chr", "", as.character(dmr$seqnames[i])) &
            RSanno$pos >= as.numeric(dmr$start[i]) &
            RSanno$pos <= as.numeric(dmr$end[i])
        out <- c(out, rownames(RSanno)[which(sel)])
    }
    unique(out)
}