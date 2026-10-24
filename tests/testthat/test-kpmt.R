test_that("KPMT p-values are between 0 and 1 and monotone", {
    set.seed(11)
    pop <- seq_len(2000)
    names(pop) <- paste0("g", seq_len(2000))
    low <- paste0("g", 1:50)
    high <- paste0("g", 1951:2000)
    p_low <- ChAMP::champ_kpmt(pop, low, tail = "lower")
    p_high <- ChAMP::champ_kpmt(pop, high, tail = "lower")
    expect_true(p_low$p.value >= 0 && p_low$p.value <= 1)
    expect_lt(p_low$p.value, 0.01)
    expect_gt(p_high$p.value, 0.9)
    expect_equal(ncol(p_low), 7L)
    expect_true("p.value" %in% colnames(p_low))
})

test_that("KPMT agrees with an exact enumeration for a tiny population", {
    N <- 14
    s <- 4
    pop <- seq_len(N)
    names(pop) <- paste0("g", seq_len(N))
    ## exhaustive enumeration of all subsets of size s
    cmb <- utils::combn(N, s)
    med <- apply(cmb, 2, function(idx) sort(idx)[ceiling(s / 2)])
    for (k in seq_len(3)) {
        obs <- paste0("g", cmb[, 10 + k])
        w <- sort(cmb[, 10 + k])[ceiling(s / 2)]
        expected <- mean(med <= w)
        got <- ChAMP::champ_kpmt(pop, obs, tail = "lower")$p.value
        expect_equal(got, expected, tolerance = 1e-8)
    }
})

test_that("KPMT rejects observations outside the population", {
    pop <- c(a = 1, b = 2, c = 3, d = 4, e = 5, f = 6)
    expect_error(ChAMP::champ_kpmt(pop, c("a", "zz")), "not found")
})
