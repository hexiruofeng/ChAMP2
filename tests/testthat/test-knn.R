test_that("KNN imputation recovers missing values", {
    set.seed(1)
    x <- matrix(stats::rnorm(200), 20, 10,
                dimnames = list(paste0("p", 1:20), paste0("s", 1:10)))
    x[sample(length(x), 20)] <- NA
    y <- ChAMP::champ_knn_impute(x, k = 3, use_impute = FALSE, strict = FALSE)
    expect_false(anyNA(y))
    expect_equal(dim(y), dim(x))
    expect_equal(dimnames(y), dimnames(x))
})

test_that("KNN imputation leaves complete matrices untouched", {
    x <- matrix(1:20, 5, 4, dimnames = list(letters[1:5], LETTERS[1:4]))
    expect_identical(ChAMP::champ_knn_impute(x, use_impute = FALSE), x)
})

test_that("rows with too many missing values get column means", {
    x <- matrix(0, 5, 4, dimnames = list(letters[1:5], LETTERS[1:4]))
    x[1, ] <- NA
    x[2, 1] <- NA
    y <- suppressMessages(ChAMP::champ_knn_impute(x, k = 2, rowmax = 0.5,
                                                  use_impute = FALSE, strict = FALSE))
    expect_false(anyNA(y[1, ]))
})

test_that("the impute package, when present, gives the same interface", {
    skip_if_not_installed("impute")
    set.seed(2)
    x <- matrix(stats::rbeta(400, 2, 5), 40, 10,
                dimnames = list(paste0("p", 1:40), paste0("s", 1:10)))
    x[cbind(sample(40, 10), sample(10, 10))] <- NA
    y1 <- ChAMP::champ_knn_impute(x, k = 5, use_impute = TRUE)
    y2 <- ChAMP::champ_knn_impute(x, k = 5, use_impute = FALSE)
    expect_false(anyNA(y1))
    expect_false(anyNA(y2))
})
