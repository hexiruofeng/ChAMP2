test_that("champ.filter removes failing probes and samples", {
    skip_if_no_champdata()
    set.seed(9)
    pf <- ChAMP:::.champ_probe_features("450K")
    ids <- rownames(pf)[seq_len(300)]
    beta <- matrix(stats::rbeta(300 * 6, 2, 5), 300, 6,
                   dimnames = list(ids, paste0("S", 1:6)))
    detP <- matrix(0.001, 300, 6, dimnames = dimnames(beta))
    detP[1:10, ] <- 0.5
    res <- suppressMessages(ChAMP::champ.filter(
        beta = beta, detP = detP, autoimpute = FALSE, filterBeads = FALSE,
        filterNoCG = FALSE, filterSNPs = TRUE, filterMultiHit = TRUE,
        filterXY = TRUE, arraytype = "450K"))
    expect_lt(nrow(res$beta), 300)
    expect_false(any(is.na(res$beta)))
    expect_true(all(res$beta > 0 & res$beta < 1))
})

test_that("champ.filter validates its inputs", {
    expect_error(suppressMessages(ChAMP::champ.filter()), "At least one")
    b <- matrix(0.5, 4, 2, dimnames = list(paste0("cg", 1:4), c("a", "b")))
    m <- b[-1, ]
    expect_error(suppressMessages(ChAMP::champ.filter(beta = b, M = m)),
                 "different rownames")
})
