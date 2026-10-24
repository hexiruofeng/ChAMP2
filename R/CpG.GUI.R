# =============================================================================
#  CpG.GUI -- interactive annotation summary for a set of CpGs
#
#  Shows the chromosome, CpG-island, feature and probe-type composition of the
#  supplied CpG list, plus the annotation table itself.
# =============================================================================

utils::globalVariables(c("Class", "Number", "Proportion"))

#' Explore the annotation of a CpG list
#'
#' @param CpG character vector of CpG identifiers. NULL uses every probe of the
#'   array.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param verbose print progress messages.
#' @param runApp launch the application; FALSE returns the shiny.appobj.
#'
#' @return invisibly, the shiny.appobj returned by \code{shiny::runApp}.
#' @export
CpG.GUI <- function(CpG = NULL, arraytype = "450K", verbose = TRUE, runApp = TRUE) {
    talk <- function(...) if (isTRUE(verbose)) .champ_msg(...)
    talk("[=============================]")
    talk("[<<<<<< CpG.GUI START >>>>>>>]")
    .champ_gui_need(dt = TRUE)
    at <- .champ_arraytype(arraytype)
    pf <- .champ_probe_features(at)
    if (is.null(CpG)) {
        CpG <- rownames(pf)
        talk("  No CpG list supplied; using all ", length(CpG), " probes.")
    }
    CpG <- unique(as.character(CpG))
    CpG <- CpG[!is.na(CpG)]
    hit <- intersect(CpG, rownames(pf))
    if (!length(hit)) .champ_stop("None of the CpGs are present in the ", at, " annotation.")
    if (length(hit) < length(CpG)) {
        talk("  ", length(CpG) - length(hit), " CpG(s) are not on the array and are ignored.")
    }
    ann <- pf[hit, , drop = FALSE]
    talk("  Summarising ", nrow(ann), " CpGs.")

    counts <- function(values) {
        values <- as.character(values)
        values[is.na(values) | !nzchar(values)] <- "NA"
        tab <- sort(table(values), decreasing = TRUE)
        data.frame(Class = names(tab), Number = as.integer(tab),
                   Proportion = as.numeric(tab) / sum(tab), stringsAsFactors = FALSE)
    }

    bar_plot <- function(df, title) {
        p <- plotly::plot_ly(df, x = ~Class, y = ~Number, type = "bar",
                             text = ~paste0(Number, " (", round(Proportion * 100, 1), "%)"),
                             textposition = "auto")
        plotly::layout(p, title = title, xaxis = list(title = ""),
                       yaxis = list(title = "number of CpGs"),
                       margin = list(b = 90, t = 50))
    }

    ui <- shiny::fluidPage(
        theme = .champ_gui_theme(),
        shiny::tags$head(shiny::tags$style(
            "#chr, #cgi, #feature, #type {height:42vh !important;}")),
        shiny::titlePanel(paste0("CpG Overview (", nrow(ann), " probes)")),
        shiny::tabsetPanel(
            shiny::tabPanel("Chromosome", plotly::plotlyOutput("chr")),
            shiny::tabPanel("CpG island", plotly::plotlyOutput("cgi")),
            shiny::tabPanel("Feature", plotly::plotlyOutput("feature")),
            shiny::tabPanel("Probe type", plotly::plotlyOutput("type")),
            shiny::tabPanel("Annotation", DT::dataTableOutput("table"))
        )
    )

    server <- function(input, output, session) {
        output$chr <- plotly::renderPlotly(bar_plot(counts(ann$CHR), "Chromosome"))
        output$cgi <- plotly::renderPlotly(bar_plot(counts(ann$cgi), "CpG island relation"))
        output$feature <- plotly::renderPlotly(bar_plot(counts(ann$feature), "Genomic feature"))
        output$type <- plotly::renderPlotly(bar_plot(counts(ann$Type), "Infinium design"))
        output$table <- DT::renderDataTable({
            cols <- intersect(c("CHR", "MAPINFO", "gene", "feature", "cgi", "Type"),
                              colnames(ann))
            DT::datatable(data.frame(Probe = rownames(ann), ann[, cols, drop = FALSE],
                                     check.names = FALSE),
                          rownames = FALSE, options = list(pageLength = 20))
        })
    }

    app <- shiny::shinyApp(ui = ui, server = server)
    talk("[<<<<<< CpG.GUI END >>>>>>>>>]")
    if (!isTRUE(runApp)) return(app)
    invisible(shiny::runApp(app))
}
