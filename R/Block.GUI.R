# =============================================================================
#  Block.GUI -- interactive explorer for champ.Block() results
#
#  Shows the block table and, for the selected block, the average beta profile
#  of the probes it contains together with the block statistics reported by
#  bumphunter.
# =============================================================================

utils::globalVariables(c("Position", "Beta", "Group", "stat", "value"))

#' Explore differentially methylated blocks in a Shiny application
#'
#' @param Block data.frame (or the list returned by champ.Block()).
#' @param beta beta value matrix used for the block analysis.
#' @param pheno phenotype vector, one entry per column of \code{beta}.
#' @param arraytype "450K", "EPICv1" or "EPICv2".
#' @param verbose print progress messages.
#' @param runApp launch the application; FALSE returns the shiny.appobj.
#'
#' @return invisibly, the shiny.appobj returned by \code{shiny::runApp}.
#' @export
Block.GUI <- function(Block = NULL, beta = NULL, pheno = NULL, arraytype = "450K",
                      verbose = TRUE, runApp = TRUE) {
    talk <- function(...) if (isTRUE(verbose)) .champ_msg(...)
    talk("[=============================]")
    talk("[<<<<<< Block.GUI START >>>>>]")
    .champ_gui_need(dt = TRUE)
    at <- .champ_arraytype(arraytype)
    beta <- .champ_as_beta(beta)

    if (is.list(Block) && !is.data.frame(Block)) {
        if (!is.null(Block$Block)) {
            Block <- Block$Block
        } else if (!length(Block)) {
            .champ_stop("'Block' is an empty list; run champ.Block() first.")
        } else {
            Block <- Block[[1L]]
        }
    }
    Block <- as.data.frame(Block, stringsAsFactors = FALSE)
    if (!nrow(Block)) .champ_stop("'Block' has no rows; nothing to display.")
    if (!"chr" %in% colnames(Block)) {
        alt <- intersect(c("seqnames", "CHR"), colnames(Block))
        if (!length(alt)) .champ_stop("'Block' has no chromosome column.")
        colnames(Block)[colnames(Block) == alt[1L]] <- "chr"
    }
    if (!all(c("start", "end") %in% colnames(Block))) {
        .champ_stop("'Block' must contain 'start' and 'end' columns.")
    }
    Block$chr <- as.character(Block$chr)
    Block$start <- as.numeric(Block$start)
    Block$end <- as.numeric(Block$end)
    Block$label <- paste0(Block$chr, ":", Block$start, "-", Block$end)

    pheno <- .champ_gui_group(pheno, ncol(beta), 4)
    pf <- .champ_probe_features(at)
    n.gr <- nlevels(pheno)
    talk("  ", nrow(Block), " blocks, ", ncol(beta), " samples.")

    region_probes <- function(i) {
        if (length(i) != 1L || is.na(i)) return(character(0))
        chr <- .champ_norm_chr(Block$chr[i])
        sel <- .champ_norm_chr(pf$CHR) == chr & !is.na(pf$MAPINFO) &
            pf$MAPINFO >= Block$start[i] & pf$MAPINFO <= Block$end[i]
        intersect(rownames(pf)[which(sel)], rownames(beta))
    }

    profile_plot <- function(i) {
        ids <- region_probes(i)
        if (!length(ids)) return(.champ_gui_empty("No probe of this block is in 'beta'."))
        ids <- ids[order(pf[ids, "MAPINFO"])]
        mat <- beta[ids, , drop = FALSE]
        pal <- .champ_gui_palette(n.gr)
        p <- plotly::plot_ly()
        for (g in seq_len(n.gr)) {
            cols <- which(pheno == levels(pheno)[g])
            if (!length(cols)) next
            p <- plotly::add_trace(
                p, x = pf[ids, "MAPINFO"],
                y = rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE),
                type = "scatter", mode = "lines+markers", name = levels(pheno)[g],
                line = list(color = pal[g], width = 3), marker = list(size = 5),
                text = ids, hovertemplate = "%{text}<br>%{x}<br>%{y:.2f}<extra></extra>")
        }
        plotly::layout(p, title = Block$label[i],
                       xaxis = list(title = "position (bp)"),
                       yaxis = list(title = "mean beta", range = c(0, 1)),
                       margin = list(b = 60, t = 50))
    }

    stats_plot <- function(i) {
        cols <- intersect(c("value", "area", "p.valueArea", "p.value", "L", "cluster"),
                          colnames(Block))
        cols <- setdiff(cols, c("start", "end"))
        if (!length(cols)) return(.champ_gui_empty("No block statistic available."))
        d <- data.frame(stat = cols,
                        value = vapply(cols, function(nm) {
                            v <- suppressWarnings(as.numeric(Block[i, nm]))
                            if (is.na(v)) 0 else v
                        }, numeric(1)))
        p <- plotly::plot_ly(d, x = ~stat, y = ~value, type = "bar",
                             text = signif(d$value, 3), textposition = "auto")
        plotly::layout(p, title = "Block statistics", xaxis = list(title = ""))
    }

    probe_table <- function(i) {
        ids <- region_probes(i)
        if (!length(ids)) return(DT::datatable(data.frame(Message = "No probe in this block.")))
        ann <- pf[ids, intersect(c("CHR", "MAPINFO", "gene", "feature", "cgi", "Type"),
                                 colnames(pf)), drop = FALSE]
        mat <- beta[ids, , drop = FALSE]
        DT::datatable(data.frame(Probe = ids, ann,
                                 MeanBeta = round(rowMeans(mat, na.rm = TRUE), 3),
                                 check.names = FALSE),
                      rownames = FALSE, options = list(pageLength = 15))
    }

    ui <- shiny::fluidPage(
        theme = .champ_gui_theme(),
        shiny::tags$head(shiny::tags$style(
            "#blocktable {font-size:80%;}",
            "#profile {height:60vh !important;}",
            "#stats {height:40vh !important;}")),
        shiny::titlePanel("Block Overview"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                width = 3,
                shiny::selectInput("block", "Block:", choices = stats::setNames(seq_len(nrow(Block)), Block$label),
                                   selected = 1L)),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel("Blocks", DT::dataTableOutput("blocktable")),
                    shiny::tabPanel("Profile",
                                    plotly::plotlyOutput("profile"),
                                    plotly::plotlyOutput("stats")),
                    shiny::tabPanel("Probes", DT::dataTableOutput("probes")))
            )
        )
    )

    server <- function(input, output, session) {
        selected <- shiny::reactive({
            i <- suppressWarnings(as.integer(input$block))
            if (is.na(i) || i < 1L || i > nrow(Block)) 1L else i
        })
        shiny::observeEvent(input$blocktable_rows_selected, {
            i <- input$blocktable_rows_selected
            if (length(i) == 1L) shiny::updateSelectInput(session, "block", selected = i)
        }, ignoreNULL = TRUE)
        output$blocktable <- DT::renderDataTable({
            DT::datatable(Block, rownames = FALSE, selection = "single",
                          options = list(pageLength = 15))
        })
        output$profile <- plotly::renderPlotly(profile_plot(selected()))
        output$stats <- plotly::renderPlotly(stats_plot(selected()))
        output$probes <- DT::renderDataTable(probe_table(selected()))
    }

    app <- shiny::shinyApp(ui = ui, server = server)
    talk("[<<<<<< Block.GUI END >>>>>>>]")
    if (!isTRUE(runApp)) return(app)
    invisible(shiny::runApp(app))
}
