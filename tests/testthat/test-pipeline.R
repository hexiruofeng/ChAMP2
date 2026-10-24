test_that("champ.DMP finds the planted difference", {
    skip_if_no_champdata()
    skip_if_not_installed("limma")
    beta <- make_beta(500L, 8L)
    pheno <- attr(beta, "group")
    dmp <- quiet_champ(champ.DMP(beta, pheno = pheno, arraytype = "450K",
                                 verbose = FALSE))
    expect_true(length(dmp) >= 1)
    tab <- dmp[[1]]
    expect_true(all(c("logFC", "P.Value", "adj.P.Val", "CHR", "MAPINFO", "gene") %in%
                        colnames(tab)))
    expect_true(any(grepl("_AVG$", colnames(tab))))
    expect_true("deltaBeta" %in% colnames(tab))
    ## the first 50 probes were shifted upwards in group B
    expect_gt(mean(tab$logFC[seq_len(20)]), 0)
})

test_that("champ.DMP handles a numeric phenotype", {
    skip_if_no_champdata()
    beta <- make_beta(300L, 8L)
    set.seed(2)
    age <- stats::rnorm(8, 50, 10)
    ## make the first 40 probes strongly associated with age
    for (i in seq_len(40)) {
        beta[i, ] <- pmin(0.99, pmax(0.01, 0.2 + 0.01 * (age - min(age)) +
                                         stats::rnorm(length(age), 0, 0.01)))
    }
    dmp <- quiet_champ(champ.DMP(beta, pheno = age, arraytype = "450K", verbose = FALSE))
    expect_true("NumericVariable" %in% names(dmp))
})

test_that("champ.QC writes its plots", {
    skip_if_no_champdata()
    beta <- make_beta(300L, 6L)
    pheno <- attr(beta, "group")
    dir <- tempfile("qc")
    quiet_champ(champ.QC(beta, pheno = pheno, PDFplot = TRUE, Rplot = FALSE,
                         resultsDir = dir, verbose = FALSE))
    expect_true(all(c("raw_mdsPlot.pdf", "raw_densityPlot.pdf",
                      "raw_SampleCluster.pdf") %in% list.files(dir)))
})

test_that("champ.SVD returns a p-value matrix", {
    skip_if_no_champdata()
    beta <- make_beta(300L, 8L)
    pheno <- attr(beta, "group")
    pd <- data.frame(Sample_Name = colnames(beta), Sample_Group = pheno,
                     stringsAsFactors = FALSE)
    dir <- tempfile("svd")
    pv <- quiet_champ(champ.SVD(beta, pd = pd, PDFplot = FALSE, Rplot = FALSE,
                                resultsDir = dir, verbose = FALSE))
    expect_true(is.matrix(pv))
    expect_true(all(pv >= 0 & pv <= 1, na.rm = TRUE))
    expect_true("Sample_Group" %in% colnames(pv))
})

test_that("champ.norm PBC keeps beta values in range", {
    skip_if_no_champdata()
    beta <- make_beta(400L, 6L)
    out <- quiet_champ(champ.norm(beta, method = "PBC", arraytype = "450K",
                                  cores = 1, verbose = FALSE,
                                  resultsDir = tempfile("norm")))
    expect_equal(dim(out), dim(beta))
    expect_true(all(out >= 0 & out <= 1))
})

test_that("champ.impute handles all three methods", {
    beta <- make_beta(200L, 6L)
    beta[sample(length(beta), 15)] <- NA
    for (m in c("Combine", "Delete", "KNN")) {
        out <- suppressMessages(champ.impute(beta, method = m, verbose = FALSE))
        expect_false(anyNA(out))
    }
})

test_that("champ.filter is idempotent on clean data", {
    skip_if_no_champdata()
    beta <- make_beta(200L, 6L)
    detP <- matrix(1e-4, nrow(beta), ncol(beta), dimnames = dimnames(beta))
    out <- suppressMessages(champ.filter(beta = beta, detP = detP,
                                         autoimpute = FALSE, filterSNPs = FALSE,
                                         filterMultiHit = FALSE, filterXY = FALSE,
                                         filterNoCG = FALSE, filterBeads = FALSE,
                                         arraytype = "450K", verbose = FALSE))
    expect_equal(rownames(out$beta), rownames(beta))
})
test_that("champ.DMR with Bumphunter finds the planted region", {
    skip_if_no_champdata()
    skip_if_not_installed("bumphunter")
    pf <- ChAMP:::.champ_probe_features("450K")
    chr1 <- pf[pf$CHR == "1" & !is.na(pf$MAPINFO), ]
    chr1 <- chr1[order(chr1$MAPINFO), ][1:800, ]
    ids <- rownames(chr1)
    set.seed(11)
    grp <- rep(c("A", "B"), each = 5)
    beta <- matrix(stats::rbeta(length(ids) * 10, 2, 8), nrow = length(ids),
                   dimnames = list(ids, paste0("S", 1:10)))
    beta[300:380, grp == "B"] <- pmin(1 - 1e-6, beta[300:380, grp == "B"] + 0.5)
    out <- quiet_champ(champ.DMR(beta, pheno = grp, method = "Bumphunter",
                                 arraytype = "450K", minProbes = 5,
                                 B = 50, cores = 1, verbose = FALSE))
    expect_true("BumphunterDMR" %in% names(out))
    expect_true(nrow(out$BumphunterDMR) >= 1)
    expect_true(all(c("seqnames", "start", "end", "width", "strand",
                      "p.valueArea") %in% colnames(out$BumphunterDMR)))
})
