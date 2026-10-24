test_that("champ.ebGSEA returns global test results", {
    skip_if_no_champdata()
    skip_if_not_installed("globaltest")
    pf <- ChAMP:::.champ_probe_features("450K")
    ids <- rownames(pf)[!is.na(pf$MAPINFO)][1:400]
    set.seed(3)
    beta <- matrix(stats::rbeta(length(ids) * 8, 2, 5), nrow = length(ids),
                   dimnames = list(ids, paste0("S", 1:8)))
    grp <- rep(c("A", "B"), each = 4)
    beta[1:40, grp == "B"] <- pmin(0.99, beta[1:40, grp == "B"] + 0.3)
    res <- quiet_champ(champ.ebGSEA(beta, pheno = grp, verbose = FALSE))
    expect_true(all(c("GSEA", "EnrichGene", "gtResult") %in% names(res)))
    expect_true(nrow(res$gtResult) > 0)
    expect_true(is.data.frame(res$gtResult))
})

test_that("champ.GSEA explains when the array subset is too small", {
    skip_if_no_champdata()
    pf <- ChAMP:::.champ_probe_features("450K")
    ids <- rownames(pf)[!is.na(pf$MAPINFO)][1:200]
    set.seed(4)
    beta <- matrix(stats::rbeta(length(ids) * 6, 2, 5), nrow = length(ids),
                   dimnames = list(ids, paste0("S", 1:6)))
    expect_error(quiet_champ(champ.GSEA(beta, Genelist = "TP53", method = "fisher",
                                         adjPval = 1, verbose = FALSE)),
                 "60%")
})
