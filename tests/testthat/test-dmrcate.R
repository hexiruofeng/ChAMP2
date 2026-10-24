test_that("the DMRcate implementation finds a planted region", {
    skip_if_no_champdata()
    skip_if_not_installed("limma")
    pf <- ChAMP:::.champ_probe_features("450K")
    ## take one chromosome, one contiguous chunk
    chr1 <- pf[pf$CHR == "1" & !is.na(pf$MAPINFO), ]
    chr1 <- chr1[order(chr1$MAPINFO), ]
    chr1 <- chr1[seq_len(1500), ]
    ids <- rownames(chr1)
    set.seed(5)
    grp <- rep(c("A", "B"), each = 5)
    beta <- matrix(stats::rbeta(length(ids) * 10, 2, 8), nrow = length(ids),
                   dimnames = list(ids, paste0("S", 1:10)))
    beta[800:860, grp == "B"] <- pmin(1 - 1e-6, beta[800:860, grp == "B"] + 0.55)
    dmr <- suppressMessages(ChAMP::champ_dmrcate(beta, grp, arraytype = "450K",
                                                 fdr = 0.05, min.cpgs = 5,
                                                 lambda = 1000, C = 2, verbose = FALSE))
    expect_true(is.data.frame(dmr))
    expect_true(nrow(dmr) >= 1)
    expect_true(all(c("seqnames", "start", "end", "no.cpgs", "Stouffer",
                      "HMFDR", "Fisher", "meandiff") %in% colnames(dmr)))
    expect_true(any(dmr$start <= chr1$MAPINFO[830] & dmr$end >= chr1$MAPINFO[830]))
})