# =============================================================================
#  champ.ebGSEA -- empirical Bayes gene set enrichment analysis
#
#  A global test (globaltest package) is applied to every gene, genes are
#  ranked by their p-value and gene sets are then tested against that ranking
#  with a Wilcoxon test and the known-population median test.
# =============================================================================

#' Gene sets from the probe annotation
#' @keywords internal
.champ_map_eid_to_cpg <- function(beta, arraytype = "450K", use_alias = TRUE) {
    pf <- .champ_probe_features(arraytype)
    RSanno <- pf[rownames(beta), c("CHR", "MAPINFO", "gene", "feature"), drop = FALSE]
    colnames(RSanno) <- c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group")
    RSanno$UCSC_RefGene_Name[is.na(RSanno$UCSC_RefGene_Name)] <- ""
    RSanno$UCSC_RefGene_Group[is.na(RSanno$UCSC_RefGene_Group)] <- ""

    ann.keep <- RSanno[substr(rownames(RSanno), 1L, 2L) == "cg" &
                           RSanno$UCSC_RefGene_Name != "", , drop = FALSE]
    geneslist <- strsplit(ann.keep$UCSC_RefGene_Name, ";", fixed = TRUE)
    names(geneslist) <- rownames(ann.keep)
    grouplist <- strsplit(ann.keep$UCSC_RefGene_Group, ";", fixed = TRUE)
    names(grouplist) <- rownames(ann.keep)

    flat <- data.frame(symbol = unlist(geneslist, use.names = FALSE),
                       group = unlist(grouplist, use.names = FALSE),
                       cpg = rep(names(geneslist), lengths(geneslist)),
                       stringsAsFactors = FALSE)
    flat$symbol <- trimws(flat$symbol)
    flat <- flat[nzchar(flat$symbol), , drop = FALSE]

    if (isTRUE(use_alias) && .champ_has("org.Hs.eg.db") && .champ_has("limma")) {
        flat$alias <- .champ_try(limma::alias2SymbolTable(flat$symbol), default = flat$symbol)
        if (length(flat$alias) != nrow(flat)) flat$alias <- flat$symbol
    } else {
        flat$alias <- flat$symbol
    }
    flat$alias[is.na(flat$alias)] <- flat$symbol[is.na(flat$alias)]
    flat <- flat[!duplicated(paste(flat$cpg, flat$alias, sep = ".")), , drop = FALSE]
    split(flat$cpg, flat$alias)
}

#' Empirical Bayes gene set enrichment
#'
#' @param beta beta value matrix.
#' @param pheno phenotype vector (binary factor or numeric covariate).
#' @param minN minimum number of genes per set.
#' @param adjPval adjusted p-value cutoff.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param cores number of workers.
#' @param verbose print progress.
#'
#' @return list with GSEA, EnrichGene and gtResult.
#' @export
champ.ebGSEA <- function(beta = NULL, pheno = NULL, minN = 5L, adjPval = 0.05,
                         arraytype = "450K", cores = 1L, verbose = TRUE) {
    .champ_require("globaltest", "champ.ebGSEA()")
    .champ_require("limma", "champ.ebGSEA()")
    beta <- .champ_as_beta(beta)
    at <- .champ_arraytype(arraytype)
    pheno <- .champ_pheno(pheno, ncol(beta))
    if (anyNA(beta)) {
        .champ_msg("champ.ebGSEA() requires complete data; imputing NA values.")
        beta <- champ_knn_impute(beta, k = 5L, strict = FALSE, verbose = FALSE)
    }

    .champ_msg("\n[ Section 1: Generate annotation ]\n")
    mapEID <- .champ_map_eid_to_cpg(beta, at)
    .champ_msg("  ", length(mapEID), " gene sets built from the annotation.")

    .champ_msg("\n[ Section 2: Global test ]\n")
    response <- if (.champ_is_numeric(pheno)) as.numeric(pheno) else
        as.numeric(factor(pheno)) - 1
    model <- if (.champ_is_numeric(pheno)) "linear" else "logistic"
    .champ_msg("  Using a ", model, " global test; this can be slow.")
    gt.o <- .champ_try(globaltest::gt(response = response,
                                      alternative = t(beta),
                                      model = model,
                                      directional = FALSE,
                                      standardize = FALSE,
                                      permutations = 0,
                                      subsets = mapEID,
                                      trace = FALSE), default = NULL)
    if (is.null(gt.o)) {
        .champ_stop("globaltest::gt() failed. Check that pheno has two groups ",
                    "(or is numeric) and that beta has no constant probes.")
    }
    resGT.m <- as.data.frame(globaltest::result(gt.o))
    if (!ncol(resGT.m)) .champ_stop("The global test returned an empty result.")
    ord <- order(resGT.m[, 1L])
    sresGT.m <- resGT.m[ord, , drop = FALSE]

    .champ_msg("\n[ Section 3: GSEA on pathways ]\n")
    PathwayList <- .champ_pathway_list()

    gseaWTfn <- function(termEID.v, rankEID.v, minN) {
        commonEID.v <- intersect(termEID.v, rankEID.v)
        nrep <- length(commonEID.v)
        if (nrep < minN) return(c(nrep, 0, 1, 1))
        otherEID.v <- setdiff(rankEID.v, termEID.v)
        if (!length(otherEID.v)) return(c(nrep, 0, 1, 1))
        rank1.idx <- match(commonEID.v, rankEID.v)
        rank2.idx <- match(otherEID.v, rankEID.v)
        w <- .champ_try(stats::wilcox.test(rank1.idx, rank2.idx, alternative = "less"),
                        default = NULL)
        if (is.null(w)) return(c(nrep, 0, 1, 1))
        n1n2 <- nrep * length(otherEID.v)
        auc <- (n1n2 - unname(w$statistic)) / n1n2
        pop.v <- seq_along(rankEID.v)
        names(pop.v) <- rankEID.v
        pvKPMT <- .champ_try(champ_kpmt(pop = pop.v, obs = commonEID.v,
                                        tail = "lower")$p.value, default = NA_real_)
        c(nrep, auc, w$p.value, pvKPMT)
    }

    rankEID.v <- rownames(sresGT.m)
    cores <- .champ_cores(cores)
    if (cores > 1L) {
        .champ_msg("  ", cores, " cores will be used.")
        cl <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        res <- parallel::parLapply(cl, PathwayList, function(x) {
            gseaWTfn(x, rankEID.v, minN)
        })
        gseaWT.m <- do.call(rbind, res)
    } else {
        gseaWT.m <- do.call(rbind, lapply(PathwayList, function(x)
            gseaWTfn(x, rankEID.v, minN)))
    }
    colnames(gseaWT.m) <- c("nREP", "AUC", "P", "P(KPMT)")
    rownames(gseaWT.m) <- names(PathwayList)

    .champ_msg("  Adjusting pathway p-values (BH).")
    ord <- order(gseaWT.m[, 3L])
    sgseaWT.m <- gseaWT.m[ord, , drop = FALSE]
    padj.v <- stats::p.adjust(sgseaWT.m[, 3L], method = "BH")
    sel.idx <- which(padj.v <= adjPval)
    .champ_msg("  ", length(sel.idx), " pathways below the adjusted p-value ", adjPval, ".")

    if (length(sel.idx) > 1L) {
        top <- cbind(sgseaWT.m[sel.idx, , drop = FALSE], padj.v[sel.idx])
        colnames(top) <- c("nREP", "AUC", "P(WT)", "P(KPMT)", "adjP")
        out <- list(top,
                    sgseaWT.m[order(sgseaWT.m[, 2L], decreasing = TRUE), , drop = FALSE])
        names(out) <- c("Rank(P)", "Rank(AUC)")
    } else if (length(sel.idx) == 1L) {
        top <- as.vector(c(sgseaWT.m[sel.idx, ], padj.v[sel.idx]))
        names(top) <- c("nREP", "AUC", "P(WT)", "P(KPMT)", "adjP")
        out <- list(top, top, rownames(sgseaWT.m)[sel.idx])
        names(out) <- c("Rank(P)", "Rank(AUC)", "POI")
    } else {
        .champ_msg("No significant pathway enriched; try relaxing adjPval.")
        out <- list()
    }
    commonEID <- lapply(PathwayList, function(x) intersect(x, rownames(sresGT.m)))
    list(GSEA = out, EnrichGene = commonEID, gtResult = sresGT.m)
}
