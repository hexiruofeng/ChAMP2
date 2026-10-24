# =============================================================================
#  QC.GUI -- interactive quality control dashboard
#
#  MDS/PCA of the samples, per-sample beta density, variance explained by the
#  principal components and the sample sheet, all rendered with plotly so the
#  plots can be zoomed and exported.
# =============================================================================

utils::globalVariables(c("Dim1", "Dim2", "Sample", "Group", "Beta", "Density",
                         "Component", "PVE"))

#' Interactive quality control
#'
#' @param beta beta value matrix (probes x samples).
#' @param pheno phenotype vector, one entry per column of \code{beta}.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param verbose print progress messages.
#' @param runApp launch the application; FALSE returns the shiny.appobj.
#' @param numPositions number of most variable probes used for the MDS/PCA.
#'
#' @return invisibly, the shiny.appobj returned by \code{shiny::runApp}.
#' @export
QC.GUI <- function(beta = NULL, pheno = NULL, arraytype = "450K", verbose = TRUE,
                   runApp = TRUE, numPositions = 1000L) {
    talk <- function(...) if (isTRUE(verbose)) .champ_msg(...)
    talk("[=============================]")
    talk("[<<<<<< QC.GUI START >>>>>>>>]")
    .champ_gui_need(dt = TRUE)
    at <- .champ_arraytype(arraytype)
    beta <- .champ_as_beta(beta)
    if (anyNA(beta)) {
        talk("  ", sum(is.na(beta)), " NA values; imputing for the QC plots.")
        beta <- champ_knn_impute(beta, k = 5L, strict = FALSE, verbose = FALSE)
    }
    pheno <- .champ_gui_group(pheno, ncol(beta), 4)
    np <- max(50L, min(as.integer(numPositions)[1L], nrow(beta)))

    v <- .champ_row_vars(beta)
    v[!is.finite(v)] <- -Inf
    idx <- order(v, decreasing = TRUE)[seq_len(np)]
    top <- beta[idx, , drop = FALSE]
    top <- top - rowMeans(top, na.rm = TRUE)
    pca <- stats::prcomp(t(top), center = FALSE, scale. = FALSE)
    pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    colnames(scores)[1:2] <- c("Dim1", "Dim2")
    scores$Sample <- rownames(scores)
    scores$Group <- pheno
    talk("  PCA on the ", np, " most variable probes.")

    mds_plot <- function() {
        pal <- .champ_gui_palette(nlevels(pheno))
        p <- plotly::plot_ly(scores, x = ~Dim1, y = ~Dim2, color = ~Group,
                             colors = pal, type = "scatter", mode = "markers",
                             text = ~Sample, marker = list(size = 9, opacity = 0.8),
                             hovertemplate = "%{text}<br>%{x:.2f}, %{y:.2f}<extra></extra>")
        plotly::layout(p, title = paste0("MDS/PCA (", np, " variable probes)"),
                       xaxis = list(title = paste0("PC1 (", round(pve[1], 1), "%)")),
                       yaxis = list(title = paste0("PC2 (", round(pve[2], 1), "%)")))
    }

    density_plot <- function() {
        n <- ncol(beta)
        if (n > 40L) {
            talk("  More than 40 samples; densities are averaged per group.")
            d <- lapply(levels(pheno), function(g) {
                cols <- which(pheno == g)
                if (!length(cols)) return(NULL)
                dd <- stats::density(rowMeans(beta[, cols, drop = FALSE], na.rm = TRUE),
                                     na.rm = TRUE)
                data.frame(Beta = dd$x, Density = dd$y, Group = g, stringsAsFactors = FALSE)
            })
            d <- do.call(rbind, d[!vapply(d, is.null, logical(1L))])
            p <- plotly::plot_ly(d, x = ~Beta, y = ~Density, color = ~Group,
                                 type = "scatter", mode = "lines")
        } else {
            d <- lapply(seq_len(n), function(j) {
                dd <- stats::density(beta[, j], na.rm = TRUE)
                data.frame(Beta = dd$x, Density = dd$y,
                           Sample = colnames(beta)[j], Group = pheno[j],
                           stringsAsFactors = FALSE)
            })
            d <- do.call(rbind, d)
            p <- plotly::plot_ly(d, x = ~Beta, y = ~Density, color = ~Sample,
                                 type = "scatter", mode = "lines")
        }
        plotly::layout(p, title = "Beta value density", xaxis = list(title = "beta"))
    }

    pve_plot <- function() {
        k <- min(20L, length(pve))
        d <- data.frame(Component = paste0("PC", seq_len(k)), PVE = pve[seq_len(k)])
        p <- plotly::plot_ly(d, x = ~Component, y = ~PVE, type = "bar")
        plotly::layout(p, title = "Variance explained", yaxis = list(title = "% variance"))
    }

    sample_table <- function() {
        d <- data.frame(Sample = colnames(beta), Group = as.character(pheno),
                        MeanBeta = round(colMeans(beta, na.rm = TRUE), 4),
                        SDBeta = round(apply(beta, 2L, stats::sd, na.rm = TRUE), 4),
                        NAs = colSums(is.na(beta)), stringsAsFactors = FALSE)
        DT::datatable(d, rownames = FALSE, options = list(pageLength = 20))
    }

    ui <- shiny::fluidPage(
        theme = .champ_gui_theme(),
        shiny::tags$head(shiny::tags$style(
            "#mds {height:55vh !important;}",
            "#density {height:55vh !important;}",
            "#pve {height:40vh !important;}")),
        shiny::titlePanel("QC Overview"),
        shiny::tabsetPanel(
            shiny::tabPanel("MDS / PCA", plotly::plotlyOutput("mds")),
            shiny::tabPanel("Density", plotly::plotlyOutput("density")),
            shiny::tabPanel("Variance", plotly::plotlyOutput("pve")),
            shiny::tabPanel("Samples", DT::dataTableOutput("samples"))
        )
    )

    server <- function(input, output, session) {
        output$mds <- plotly::renderPlotly(mds_plot())
        output$density <- plotly::renderPlotly(density_plot())
        output$pve <- plotly::renderPlotly(pve_plot())
        output$samples <- DT::renderDataTable(sample_table())
    }

    app <- shiny::shinyApp(ui = ui, server = server)
    talk("[<<<<<< QC.GUI END >>>>>>>>]")
    if (!isTRUE(runApp)) return(app)
    invisible(shiny::runApp(app))
}