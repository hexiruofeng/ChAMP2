# =============================================================================
#  champ.DMP -- differentially methylated positions
# =============================================================================

#' Differential methylation analysis at the probe level
#'
#' For a categorical phenotype every pair of groups is compared with limma
#' (moderated t-test); for a numeric phenotype a linear regression is fitted.
#' Results are annotated with the ChAMP probe annotation.
#'
#' @param beta beta value matrix (probes x samples).
#' @param pheno phenotype vector (character, factor or numeric covariate).
#' @param compare.group two group labels to compare (optional).
#' @param adjPVal adjusted p-value cutoff used for reporting.
#' @param adjust.method multiple testing correction ("BH" by default).
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param verbose print progress.
#'
#' @return named list of data.frames, one per comparison (or "NumericVariable").
#' @export
champ.DMP <- function(beta = NULL, pheno = NULL, compare.group = NULL,
                      adjPVal = 0.05, adjust.method = "BH",
                      arraytype = "450K", verbose = TRUE) {
    .champ_msg("[===========================]")
    .champ_msg("[<<<<< ChAMP.DMP START >>>>>]")
    .champ_msg("-----------------------------")
    .champ_require("limma", "champ.DMP()")
    beta <- .champ_as_beta(beta)
    at <- .champ_arraytype(arraytype)
    pheno <- .champ_pheno(pheno, ncol(beta))
    if (length(unique(pheno)) < 2L) {
        .champ_stop("'pheno' must contain at least two groups.")
    }
    if (anyNA(beta)) {
        .champ_msg(sum(is.na(beta)), " NA values detected; imputing before the test.")
        beta <- champ_knn_impute(beta, k = 5L, strict = FALSE, verbose = FALSE)
    }

    .champ_msg("!!! champ.DMP() notes !!!\n")
    .champ_msg("    (1) With more than two groups every pair is compared; use compare.group to restrict.\n")
    .champ_msg("    (2) A numeric pheno is analysed by linear regression.\n")
    .champ_msg("--------------------------------")

    numeric_pheno <- .champ_is_numeric(pheno)
    Compare <- NULL
    if (numeric_pheno) {
        .champ_msg("  Your pheno is numeric (", length(pheno), " values, range ",
                   min(pheno), "-", max(pheno), ").")
    } else {
        .champ_msg("  Your pheno is ", class(pheno), " with the following groups:")
        for (g in sort(unique(pheno))) {
            .champ_msg("    <", g, ">: ", sum(pheno == g), " samples.")
        }
        ug <- sort(unique(pheno))
        if (length(ug) == 2L) {
            if (!is.null(compare.group) && !all(compare.group %in% ug)) {
                .champ_stop("'compare.group' is not in accord with 'pheno'.")
            }
            Compare <- list(ug)
        } else {
            if (is.null(compare.group)) {
                Compare <- as.list(as.data.frame(utils::combn(ug, 2L), stringsAsFactors = FALSE))
            } else if (length(compare.group) == 2L && all(compare.group %in% ug)) {
                Compare <- list(sort(compare.group))
            } else {
                .champ_stop("'compare.group' must contain two values found in 'pheno'.")
            }
        }
        names(Compare) <- vapply(Compare, function(x) paste(x[1L], x[2L], sep = "_to_"), character(1L))
    }

    calculate <- function(group) {
        .champ_msg("  -----------------------------")
        .champ_msg("  Comparing: ", group[1L], " vs ", group[2L])
        sel <- pheno %in% group
        p <- factor(pheno[sel], levels = group)
        if (any(table(p) < 2L)) {
            .champ_warn("Group ", paste(names(which(table(p) < 2)), collapse = ","),
                        " has fewer than 2 samples; skipping this comparison.")
            return(NULL)
        }
        design <- stats::model.matrix(~0 + p)
        colnames(design) <- levels(p)
        contrast.matrix <- limma::makeContrasts(contrasts = paste(levels(p)[2L], "-",
                                                                  levels(p)[1L]),
                                                levels = colnames(design))
        fit <- limma::lmFit(beta[, sel, drop = FALSE], design)
        fit2 <- limma::contrasts.fit(fit, contrast.matrix)
        fit3 <- tryCatch(limma::eBayes(fit2), error = function(e) {
            .champ_stop("limma failed: ", conditionMessage(e),
                        "\n  This usually means that all samples in a group have identical values.")
        })
        tab <- limma::topTable(fit3, coef = 1L, number = nrow(beta),
                               adjust.method = adjust.method, p.value = 1)
        .champ_msg("  ", sum(tab$adj.P.Val <= adjPVal), " significant DMPs at adjusted p <= ",
                   adjPVal, ".")
        tab
    }

    DMPs <- list()
    if (numeric_pheno) {
        .champ_msg("\n[ Section 2: Linear regression on a numeric covariate ]\n")
        df <- data.frame(pheno = pheno)
        design <- stats::model.matrix(~pheno, data = df)
        fit <- limma::lmFit(beta, design)
        fit <- limma::eBayes(fit)
        tab <- limma::topTable(fit, coef = 2L, number = nrow(beta),
                               adjust.method = adjust.method, p.value = 1)
        .champ_msg("  ", sum(tab$adj.P.Val <= adjPVal), " significant CpGs.\n")
        if (any(tab$adj.P.Val <= adjPVal)) DMPs[["NumericVariable"]] <- tab
    } else {
        .champ_msg("\n[ Section 2: Pairwise differential methylation ]\n")
        for (nm in names(Compare)) {
            tab <- calculate(Compare[[nm]])
            if (!is.null(tab) && any(tab$adj.P.Val <= adjPVal)) DMPs[[nm]] <- tab
        }
    }
    if (!length(DMPs)) {
        .champ_stop("champ.DMP() did not detect a single significant CpG. ",
                    "You may relax adjPVal or use a different adjustment method.")
    }

    .champ_msg("\n[ Section 3: Add annotation ]\n")
    pf <- .champ_probe_features(at)
    for (nm in names(DMPs)) {
        tab <- DMPs[[nm]]
        com.idx <- intersect(rownames(tab), rownames(pf))
        if (!length(com.idx)) next
        tab <- tab[com.idx, , drop = FALSE]
        if (!numeric_pheno && nm %in% names(Compare)) {
            g <- Compare[[nm]]
            avg <- cbind(
                rowMeans(beta[com.idx, pheno == g[1L], drop = FALSE]),
                rowMeans(beta[com.idx, pheno == g[2L], drop = FALSE]))
            colnames(avg) <- paste(g, "AVG", sep = "_")
            avg <- cbind(avg, deltaBeta = avg[, 2L] - avg[, 1L])
            DMPs[[nm]] <- data.frame(tab, avg, pf[com.idx, , drop = FALSE],
                                     check.names = FALSE)
        } else {
            DMPs[[nm]] <- data.frame(tab, pf[com.idx, , drop = FALSE], check.names = FALSE)
        }
    }
    .champ_msg("[<<<<<< ChAMP.DMP END >>>>>>]")
    .champ_msg("[===========================]")
    .champ_msg("[You may want to process DMP.GUI() or champ.GSEA() next.]\n")
    DMPs
}