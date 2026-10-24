test_that("GUIs build a Shiny app object without launching it", {
    skip_if_no_champdata()
    skip_if_not_installed("shiny")
    skip_if_not_installed("plotly")
    skip_if_not_installed("DT")
    skip_if_not_installed("shinythemes")

    beta <- make_beta(400L, 8L)
    pheno <- attr(beta, "group")
    quiet_champ({
        dmp <- champ.DMP(beta, pheno = pheno, arraytype = "450K", verbose = FALSE)
    })

    app <- suppressMessages(DMP.GUI(dmp[[1]], beta = beta, pheno = pheno,
                                    runApp = FALSE, verbose = FALSE))
    expect_s3_class(app, "shiny.appobj")

    app <- suppressMessages(QC.GUI(beta, pheno = pheno, runApp = FALSE, verbose = FALSE))
    expect_s3_class(app, "shiny.appobj")

    app <- suppressMessages(CpG.GUI(rownames(beta)[1:100], runApp = FALSE, verbose = FALSE))
    expect_s3_class(app, "shiny.appobj")
})

test_that("GUI input validation is actionable", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("plotly")
    skip_if_not_installed("DT")
    expect_error(suppressMessages(DMP.GUI(NULL)), "beta")
    expect_error(suppressMessages(CpG.GUI("not-a-probe")), "annotation")
})