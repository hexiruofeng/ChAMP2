test_that("array type normalisation accepts every documented spelling", {
    clear_champ_warnings()
    expect_equal(ChAMP:::.champ_arraytype("450K"), "450K")
    expect_equal(ChAMP:::.champ_arraytype("450k"), "450K")
    expect_equal(ChAMP:::.champ_arraytype("HM450"), "450K")
    expect_equal(ChAMP:::.champ_arraytype("EPICv2"), "EPICv2")
    expect_equal(ChAMP:::.champ_arraytype("EPICv1"), "EPICv1")
    expect_equal(ChAMP:::.champ_arraytype("EPIC-1"), "EPICv1")
    expect_equal(ChAMP:::.champ_arraytype("27K"), "27K")
    expect_error(ChAMP:::.champ_arraytype("Nope"), "Unsupported")
    expect_warning(ChAMP:::.champ_arraytype("EPIC"), "EPIC version 2")
})

test_that("genome mapping matches the annotation builds", {
    expect_equal(ChAMP:::.champ_genome("450K"), "hg19")
    expect_equal(ChAMP:::.champ_genome("EPICv1"), "hg19")
    expect_equal(ChAMP:::.champ_genome("EPICv2"), "hg38")
})

test_that("beta / M transformations round trip", {
    b <- matrix(c(0.1, 0.5, 0.9, 0.25), 2)
    m <- ChAMP:::.champ_beta2m(b)
    expect_equal(ChAMP:::.champ_m2beta(m), b, tolerance = 1e-12)
    expect_equal(ChAMP:::.champ_ilogit2(ChAMP:::.champ_logit2(b)), b, tolerance = 1e-12)
})

test_that("pheno coercion is strict about lengths", {
    expect_equal(ChAMP:::.champ_pheno(factor(c("a", "b"))), c("a", "b"))
    expect_error(ChAMP:::.champ_pheno(c("a", "b"), n = 3), "length")
})

test_that("probe annotation is normalised", {
    skip_if_no_champdata()
    pf <- ChAMP:::.champ_probe_features("450K")
    expect_true(all(c("CHR", "MAPINFO", "gene", "feature", "cgi", "Type") %in% colnames(pf)))
    expect_false(any(grepl("^chr", pf$CHR[!is.na(pf$CHR)])))
    expect_true(all(pf$cgi[!is.na(pf$cgi)] %in% c("island", "opensea", "shelf", "shore")))
})

test_that("dependency report covers every exported helper", {
    tab <- champ_check_env(quiet = TRUE)
    expect_true(all(c("package", "role", "installed") %in% colnames(tab)))
    expect_true("limma" %in% tab$package)
    expect_true("kpmt" %in% tab$package)
})