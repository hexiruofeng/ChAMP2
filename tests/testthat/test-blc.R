test_that("beta mixture EM recovers a known mixture", {
    set.seed(7)
    x <- c(stats::rbeta(400, 2, 20), stats::rbeta(300, 20, 2), stats::rbeta(300, 5, 5))
    w <- matrix(0, length(x), 3)
    w[x < 0.2, 1] <- 1
    w[x > 0.8, 3] <- 1
    w[w[, 1] == 0 & w[, 3] == 0, 2] <- 1
    fit <- ChAMP::champ_blc(matrix(x, ncol = 1), w, maxiter = 30, tol = 1e-6)
    expect_equal(dim(fit$a), c(3L, 1L))
    expect_equal(dim(fit$b), c(3L, 1L))
    expect_equal(length(fit$eta), 3L)
    expect_equal(sum(fit$eta), 1, tolerance = 1e-8)
    expect_true(all(fit$a > 0))
    expect_true(all(fit$b > 0))
    ## the low component must be centred well below the high component
    expect_lt(fit$mu[1, 1], fit$mu[3, 1])
    expect_equal(rowSums(fit$w), rep(1, length(x)), tolerance = 1e-8)
})

test_that("BMIQ produces a normalised profile in [0, 1]", {
    skip_if_no_champdata()
    pf <- ChAMP:::.champ_probe_features("450K")
    set.seed(3)
    ids <- sample(rownames(pf), 2000)
    beta <- stats::rbeta(length(ids), 2, 5)
    names(beta) <- ids
    design <- ifelse(pf[ids, "Type"] == "I", 1L, 2L)
    res <- suppressMessages(ChAMP::champ.BMIQ(beta, design, niter = 5, plots = FALSE,
                                              verbose = FALSE))
    expect_equal(length(res$nbeta), length(beta))
    expect_true(all(res$nbeta >= 0 & res$nbeta <= 1))
    expect_equal(res$nbeta[design == 1], beta[design == 1])
})
