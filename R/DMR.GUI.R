# =============================================================================
#  DMR.GUI -- interactive explorer for champ.DMR() results
#
#  Shows the DMR table, the probe-level beta profile of the selected region
#  (one line per phenotype group, one point per probe), the region statistics
#  and the annotation of the probes it contains.
# =============================================================================

utils::globalVariables(c("Position", "Beta", "Group", "Probe", "value", "stat"))

#' Normalise a DMR table coming from champ.DMR()/DMRcate
#' @keywords internal
.champ_gui_dmr_table <- function(DMR, verbose = TRUE) {
    if (is.list(DMR) && !is.data.frame(DMR)) {
        if (!length(DMR)) .champ_stop("'DMR' is an empty list; run champ.DMR() first.")
        if (isTRUE(verbose)) .champ_msg("  'DMR' is a list; using '", names(DMR)[1L], "'.")
        DMR <- DMR[[1L]]
    }
    DMR <- as.data.frame(DMR, stringsAsFactors = FALSE)
    if (!nrow(DMR)) .champ_stop("'DMR' has no rows; nothing to display.")
    if (!all(c("start", "end") %in% colnames(DMR))) {
        .champ_stop("'DMR' must contain 'start' and 'end' columns (champ.DMR() output).")
    }
    if (!"seqnames" %in% colnames(DMR)) {
        alt <- intersect(c("chr", "CHR", "chrom"), colnames(DMR))
        if (!length(alt)) .champ_stop("'DMR' has no chromosome column.")
        colnames(DMR)[colnames(DMR) == alt[1L]] <- "seqnames"
    }
    DMR$seqnames <- .champ_norm_chr(DMR$seqnames)
    DMR$start <- as.numeric(DMR$start)
    DMR$end <- as.numeric(DMR$end)
    keep <- !is.na(DMR$seqnames) & !is.na(DMR$start) & !is.na(DMR$end)
    if (!all(keep)) {
        if (isTRUE(verbose)) .champ_msg("  Dropping ", sum(!keep), " DMR(s) without coordinates.")
        DMR <- DMR[keep, , drop = FALSE]
    }
    if (!nrow(DMR)) .champ_stop("No DMR with valid coordinates.")
    DMR$label <- paste0(DMR$seqnames, ":", DMR$start, "-", DMR$end)
    DMR
}

#' Explore differentially methylated regions in a Shiny application
#'
#' @param DMR data.frame (or list) returned by champ.DMR().
#' @param beta beta value matrix used for the DMR analysis.
#' @param pheno phenotype vector, one entry per column of \code{beta}.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param verbose print progress messages.
#' @param runApp launch the application; FALSE returns the shiny.appobj.
#'
#' @return invisibly, the shiny.appobj returned by \code{shiny::runApp}.
#' @export
DMR.GUI <- function(DMR = NULL, beta = NULL, pheno = NULL, arraytype = "450K",
                    verbose = TRUE, runApp = TRUE) {
    talk <- function(...) if (isTRUE(verbose)) .champ_msg(...)
    talk("[=============================]")
    talk("[<<<<<< DMR.GUI START >>>>>>>]")
    .champ_gui_need(dt = TRUE)
    at <- .champ_arraytype(arraytype)
    beta <- .champ_as_beta(beta)
    DMR <- .champ_gui_dmr_table(DMR, verbose)
    pheno <- .champ_gui_group(pheno, ncol(beta), 4)
    pf <- .champ_probe_features(at)
    n.gr <- nlevels(pheno)
    talk("  ", nrow(DMR), " DMRs, ", ncol(beta), " samples, ", n.gr, " group(s).")

    region_probes <- function(i) {
        if (length(i) != 1L || is.na(i)) return(character(0))
        sel <- pf$CHR == DMR$seqnames[i] & !is.na(pf$MAPINFO) &
            pf$MAPINFO >= DMR$start[i] & pf$MAPINFO <= DMR$end[i]
        ids <- rownames(pf)[which(sel)]
        intersect(ids, rownames(beta))
    }

    profile_plot <- function(i) {
        ids <- region_probes(i)
        if (!length(ids)) return(.champ_gui_empty("No probe of this region is in 'beta'."))
        ids <- ids[order(pf[ids, "MAPINFO"])]
        mat <- beta[ids, , drop = FALSE]
        pal <- .champ_gui_palette(n.gr)
        p <- plotly::plot_ly()
        for (g in seq_len(n.gr)) {
            cols <- which(pheno == levels(pheno)[g])
            if (!length(cols)) next
            p <- plotly::add_trace(
                p, x = pf[ids, "MAPINFO"], y = rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE),
                type = "scatter", mode = "lines+markers", name = levels(pheno)[g],
                line = list(color = pal[g], width = 3), marker = list(size = 5),
                text = ids, hovertemplate = "%{text}<br>%{x}<br>%{y:.2f}<extra></extra>")
        }
        plotly::layout(p, title = DMR$label[i],
                       xaxis = list(title = "position (bp)"),
                       yaxis = list(title = "mean beta", range = c(0, 1)),
                       margin = list(b = 60, t = 50))
    }

    stats_plot <- function(i) {
        cols <- intersect(c("Stouffer", "HMFDR", "Fisher", "min_smoothed_fdr"),
                          colnames(DMR))
        if (!length(cols)) return(.champ_gui_empty("No region statistics in this DMR table."))
        vals <- vapply(cols, function(nm) suppressWarnings(as.numeric(DMR[i, nm])), numeric(1))
        d <- data.frame(stat = cols, value = -log10(pmax(vals, 1e-300)))
        p <- plotly::plot_ly(d, x = ~stat, y = ~value, type = "bar",
                             text = signif(vals, 3), textposition = "auto")
        plotly::layout(p, title = "Region statistics (-log10)",
                       xaxis = list(title = ""), yaxis = list(title = "-log10 p"))
    }

    probe_table <- function(i) {
        ids <- region_probes(i)
        if (!length(ids)) return(DT::datatable(data.frame(Message = "No probe in this region.")))
        ann <- pf[ids, intersect(c("CHR", "MAPINFO", "gene", "feature", "cgi", "Type"),
                                 colnames(pf)), drop = FALSE]
        DT::datatable(data.frame(Probe = ids, ann, check.names = FALSE),
                      rownames = FALSE, options = list(pageLength = 15))
    }

    ui <- shiny::fluidPage(
        theme = .champ_gui_theme(),
        shiny::tags$head(shiny::tags$style(
            "#dmrtable {font-size:80%;}",
            "#profile {height:60vh !important;}",
            "#stats {height:40vh !important;}")),
        shiny::titlePanel("DMR Overview"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                width = 3,
                shiny::selectInput("dmr", "DMR:",
                                   choices = stats::setNames(seq_len(nrow(DMR)), DMR$label),
                                   selected = 1L),
                shiny::helpText("Pick a region, or click a row in the table.")),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel("Regions",
                                    DT::dataTableOutput("dmrtable")),
                    shiny::tabPanel("Profile",
                                    plotly::plotlyOutput("profile"),
                                    plotly::plotlyOutput("stats")),
                    shiny::tabPanel("Probes",
                                    DT::dataTableOutput("probes")))
            )
        )
    )

    server <- function(input, output, session) {
        selected <- shiny::reactive({
            i <- suppressWarnings(as.integer(input$dmr))
            if (is.na(i) || i < 1L || i > nrow(DMR)) 1L else i
        })
        shiny::observeEvent(input$dmrtable_rows_selected, {
            i <- input$dmrtable_rows_selected
            if (length(i) == 1L) shiny::updateSelectInput(session, "dmr", selected = i)
        }, ignoreNULL = TRUE)
        output$dmrtable <- DT::renderDataTable({
            cols <- intersect(c("seqnames", "start", "end", "width", "no.cpgs",
                                "min_smoothed_fdr", "Stouffer", "HMFDR", "Fisher",
                                "maxdiff", "meandiff", "overlapping.genes"),
                              colnames(DMR))
            DT::datatable(DMR[, cols, drop = FALSE], rownames = FALSE,
                          selection = "single", options = list(pageLength = 15))
        })
        output$profile <- plotly::renderPlotly(profile_plot(selected()))
        output$stats <- plotly::renderPlotly(stats_plot(selected()))
        output$probes <- DT::renderDataTable(probe_table(selected()))
    }

    app <- shiny::shinyApp(ui = ui, server = server)
    talk("[<<<<<< DMR.GUI END >>>>>>>>>]")
    if (!isTRUE(runApp)) return(app)
    invisible(shiny::runApp(app))
}
