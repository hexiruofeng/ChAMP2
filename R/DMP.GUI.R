# =============================================================================
#  DMP.GUI -- interactive explorer for champ.DMP() results
#
#  Builds a Shiny application that shows the differentially methylated
#  positions of one comparison:
#    * a filterable table of the significant DMPs,
#    * a median-centred heatmap of their beta values (number of probes capped),
#    * CpG-island / feature / feature-CpG composition barplots,
#    * a gene view (probe profile with feature tracks + gene enrichment),
#    * a volcano plot built from the real logFC and adj.P.Val columns, and
#    * a single-CpG boxplot / scatter plot across the samples.
#
#  Only shiny, plotly, DT (tables) and shinythemes (theme) are needed; each is
#  loaded through .champ_require() so a missing package produces an actionable
#  error.  The .champ_gui_* helpers below are shared with DMR.GUI(),
#  Block.GUI(), QC.GUI() and CpG.GUI().
# =============================================================================

utils::globalVariables(c("pos", "Value", "Class", "Group", "Proportion", "Number",
                         "logFC", "neglog10P", "gene", "probe", "sample", "x", "y",
                         "Beta", "Pheno", "label"))

# ---------------------------------------------------------------------------
# helpers shared by all *_GUI files
# ---------------------------------------------------------------------------

.champ_gui_need <- function(dt = FALSE) {
    .champ_require("shiny", "the ChAMP graphical interfaces")
    .champ_require("plotly", "the ChAMP graphical interfaces")
    if (isTRUE(dt)) .champ_require("DT", "the interactive tables of the ChAMP GUIs")
    invisible(TRUE)
}

.champ_gui_theme <- function() {
    if (.champ_has("shinythemes")) shinythemes::shinytheme("readable") else NULL
}

.champ_gui_palette <- function(n) {
    n <- max(1L, as.integer(n)[1L])
    grDevices::hcl.colors(n, if (n <= 8L) "Dark 3" else "Viridis")
}

.champ_gui_empty <- function(msg) {
    plotly::layout(plotly::plot_ly(), title = list(text = msg, x = 0.5),
                   xaxis = list(visible = FALSE), yaxis = list(visible = FALSE))
}

#' Grouping factor used by every GUI (numeric phenotypes are binned).
.champ_gui_group <- function(pheno, n, ncut = 4L) {
    pheno <- .champ_pheno(pheno, n)
    if (length(unique(pheno[!is.na(pheno)])) < 2L) {
        .champ_stop("'pheno' must contain at least two distinct values.")
    }
    if (.champ_is_numeric(pheno)) {
        k <- max(2L, min(length(unique(pheno)), as.integer(ncut)[1L]))
        return(base::cut(pheno, breaks = k, include.lowest = TRUE, dig.lab = 3L))
    }
    factor(pheno, levels = sort(unique(pheno)))
}

#' Top-variance, median-centred matrix ready for plotly::plot_ly(type = "heatmap").
.champ_gui_heat <- function(beta, max.probes = 500L) {
    beta <- as.matrix(beta)
    if (!nrow(beta) || !ncol(beta)) return(NULL)
    v <- .champ_row_vars(beta)
    v[!is.finite(v)] <- -Inf
    n <- min(max(1L, as.integer(max.probes)[1L]), nrow(beta))
    z <- beta[order(v, decreasing = TRUE)[seq_len(n)], , drop = FALSE]
    z <- z - apply(z, 1L, stats::median, na.rm = TRUE)
    z[!is.finite(z)] <- NA_real_
    lim <- suppressWarnings(max(abs(z), na.rm = TRUE))
    if (!is.finite(lim) || lim <= 0) lim <- 0.5
    list(z = z, zlim = c(-lim, lim))
}

#' Horizontal colour tracks (feature / CpG island) drawn under a probe profile.
.champ_gui_track <- function(values, palette, y) {
    values <- as.character(values)
    values[is.na(values)] <- "NA"
    out <- list()
    for (cls in names(palette)) {
        idx <- which(values == cls)
        if (!length(idx)) next
        out[[length(out) + 1L]] <- list(
            x = as.numeric(rbind(idx - 0.5, idx + 0.5, NA_real_)),
            y = rep(c(y, y, NA_real_), length(idx)),
            name = cls, colour = unname(palette[[cls]]))
    }
    out
}

.champ_gui_add_track <- function(p, tracks) {
    for (tr in tracks) {
        p <- plotly::add_trace(p, x = tr$x, y = tr$y, type = "scatter", mode = "lines",
                               name = tr$name, text = tr$name, hoverinfo = "text",
                               line = list(color = tr$colour, width = 8), opacity = 0.6)
    }
    p
}

.champ_gui_feature_cols <- c("1stExon" = "#00eeee", "3'UTR" = "#8b008b",
                             "5'UTR" = "#00ee00", "Body" = "#ffd700", "IGR" = "#9e9e9e",
                             "TSS1500" = "#ff3030", "TSS200" = "#00688b", "NA" = "#d9d9d9")

.champ_gui_cgi_cols <- c("island" = "#ff83fa", "opensea" = "#ffdab9",
                         "shelf" = "#bbffff", "shore" = "#98fb98", "NA" = "#d9d9d9")

# ---------------------------------------------------------------------------
# DMP specific helpers
# ---------------------------------------------------------------------------

#' Validate a DMP table and fill in any missing annotation column.
.champ_gui_dmp_table <- function(DMP, arraytype, verbose = TRUE) {
    if (is.list(DMP) && !is.data.frame(DMP)) {
        if (!length(DMP)) .champ_stop("'DMP' is an empty list; run champ.DMP() first.")
        if (isTRUE(verbose)) .champ_msg("  'DMP' is a list; using '", names(DMP)[1L], "'.")
        DMP <- DMP[[1L]]
    }
    DMP <- as.data.frame(DMP, stringsAsFactors = FALSE)
    if (!nrow(DMP)) .champ_stop("'DMP' has no rows; nothing to display.")
    miss <- setdiff(c("logFC", "adj.P.Val"), colnames(DMP))
    if (length(miss)) {
        .champ_stop("'DMP' has no column(s): ", paste(miss, collapse = ", "), ".\n",
                    "  Pass DMP[[1]] from champ.DMP(), which has logFC and adj.P.Val.")
    }
    if (!.champ_is_numeric(DMP$logFC) || !.champ_is_numeric(DMP$adj.P.Val)) {
        .champ_stop("'DMP$logFC' and 'DMP$adj.P.Val' must be numeric.")
    }
    add <- setdiff(c("CHR", "MAPINFO", "gene", "feature", "cgi"), colnames(DMP))
    if (length(add)) {
        ann <- .champ_probe_features_for(rownames(DMP), arraytype)
        for (nm in add) DMP[[nm]] <- ann[[nm]]
        if (isTRUE(verbose)) {
            .champ_msg("  Added ", paste(add, collapse = ", "), " from the ChAMP annotation.")
        }
    }
    DMP$CHR <- as.character(DMP$CHR)
    DMP$feature <- as.character(DMP$feature)
    DMP$cgi <- as.character(DMP$cgi)
    DMP$gene <- as.character(DMP$gene)
    DMP$feature[is.na(DMP$feature)] <- "NA"
    DMP$cgi[is.na(DMP$cgi)] <- "NA"
    DMP$gene[is.na(DMP$gene)] <- ""
    if (is.null(DMP$feat.cgi)) DMP$feat.cgi <- paste(DMP$feature, DMP$cgi, sep = " - ")
    DMP
}

#' Explore differentially methylated positions in a Shiny application
#'
#' Opens a Shiny/plotly interface with the DMP table, heatmap, feature/CpG
#' composition, per-gene profile and per-CpG distribution of a comparison
#' produced by \code{\link{champ.DMP}}.
#'
#' @param DMP data.frame (or list of data.frames) returned by champ.DMP().
#' @param beta beta value matrix used for the DMP test (probes x samples).
#' @param pheno phenotype vector, one entry per column of \code{beta}.  A
#'   numeric phenotype is binned into \code{cutgroupnumber} groups.
#' @param cutgroupnumber number of groups used when \code{pheno} is numeric.
#' @param arraytype "450K", "EPICv1" or "EPICv2"; used for missing annotation.
#' @param verbose print progress messages.
#' @param runApp launch the application; when FALSE the shiny.appobj is
#'   returned instead (useful for tests).
#'
#' @return invisibly, the shiny.appobj returned by \code{shiny::runApp}.
#' @export
DMP.GUI <- function(DMP = NULL, beta = NULL, pheno = NULL, cutgroupnumber = 4,
                    arraytype = "450K", verbose = TRUE, runApp = TRUE) {
    talk <- function(...) if (isTRUE(verbose)) .champ_msg(...)
    talk("[=============================]")
    talk("[<<<<<< DMP.GUI START >>>>>>>]")
    .champ_gui_need(dt = TRUE)
    at <- .champ_arraytype(arraytype)
    beta <- .champ_as_beta(beta)
    DMP <- .champ_gui_dmp_table(DMP, at, verbose)
    pheno <- .champ_gui_group(pheno, ncol(beta), cutgroupnumber)
    probes <- intersect(rownames(DMP), rownames(beta))
    if (!length(probes)) {
        .champ_stop("None of the probes in 'DMP' are present in 'beta'.\n",
                    "  Make sure both objects come from the same data set.")
    }
    if (length(probes) < nrow(DMP)) {
        talk("  ", nrow(DMP) - length(probes), " DMP probe(s) are absent from 'beta'.")
    }
    DMP <- DMP[rownames(DMP) %in% probes, , drop = FALSE]
    DMP$logFC[is.na(DMP$logFC)] <- 0
    DMP$adj.P.Val[is.na(DMP$adj.P.Val)] <- 1
    n.gr <- nlevels(pheno)
    genes <- sort(table(DMP$gene[nzchar(DMP$gene)]), decreasing = TRUE)
    gene_default <- if (length(genes)) names(genes)[1L] else ""
    talk("  ", nrow(DMP), " DMPs, ", ncol(beta), " samples, ", n.gr, " phenotype group(s).")

    prop_table <- function(all_values, sel_values, lfc, fc) {
        all_values <- as.character(all_values)
        sel_values <- as.character(sel_values)
        all_values[is.na(all_values)] <- "NA"
        sel_values[is.na(sel_values)] <- "NA"
        lv <- sort(unique(c(all_values, sel_values)))
        tab <- rbind(All = table(factor(all_values, levels = lv)),
                     Hyper = table(factor(sel_values[lfc >= fc], levels = lv)),
                     Hypo = table(factor(sel_values[lfc <= -fc], levels = lv)))
        den <- rowSums(tab)
        den[den == 0] <- 1
        data.frame(Group = rep(rownames(tab), times = length(lv)),
                   Class = rep(lv, each = nrow(tab)),
                   Number = as.vector(tab),
                   Proportion = as.vector(tab) / rep(den, times = length(lv)),
                   stringsAsFactors = FALSE)
    }

    dmptable <- function(d) {
        if (!nrow(d)) return(DT::datatable(data.frame(Message = "No DMP at these cutoffs.")))
        cols <- intersect(c("CHR", "MAPINFO", "gene", "feature", "cgi", "logFC",
                            "deltaBeta", grep("_AVG$", colnames(d), value = TRUE),
                            "t", "P.Value", "adj.P.Val", "B"), colnames(d))
        DT::datatable(data.frame(ID = rownames(d), d[, cols, drop = FALSE],
                                 check.names = FALSE),
                      rownames = FALSE, options = list(pageLength = 20))
    }

    heat_plot <- function(sel) {
        if (!nrow(sel)) return(.champ_gui_empty("No significant DMP for these cutoffs."))
        h <- .champ_gui_heat(beta[rownames(sel), , drop = FALSE], max.probes = 500L)
        if (is.null(h)) return(.champ_gui_empty("No probe to display."))
        p <- plotly::plot_ly(x = colnames(h$z), y = rownames(h$z), z = h$z,
                             type = "heatmap", zmin = h$zlim[1], zmax = h$zlim[2],
                             colors = grDevices::colorRampPalette(
                                 c("deepskyblue", "black", "brown1"))(101),
                             hovertemplate = "%{y}<br>%{x}<br>%{z:.2f}<extra></extra>")
        plotly::layout(p, title = paste0("Median-centred beta of ", nrow(h$z), " probes"),
                       margin = list(l = 130, r = 30, b = 80, t = 50))
    }

    bar_plot <- function(df, title) {
        if (!nrow(df)) return(.champ_gui_empty(title))
        p <- plotly::plot_ly(df, x = ~Group, y = ~Proportion, color = ~Class,
                             type = "bar", text = ~Number,
                             hovertemplate = "%{x} %{data.name}: %{text}<extra></extra>")
        plotly::layout(p, title = title, barmode = "stack",
                       xaxis = list(title = ""), yaxis = list(title = "proportion"),
                       margin = list(l = 60, r = 20, b = 100, t = 50))
    }

    gene_enrich <- function(sel) {
        sel <- sel[nzchar(sel$gene), , drop = FALSE]
        if (!nrow(sel)) return(.champ_gui_empty("No gene annotation for these DMPs."))
        dir <- factor(ifelse(sel$logFC > 0, "Hyper", "Hypo"), levels = c("Hyper", "Hypo"))
        tab <- table(sel$gene, dir)
        keep <- utils::head(names(sort(rowSums(tab), decreasing = TRUE)), 70L)
        tab <- tab[keep, , drop = FALSE]
        df <- data.frame(gene = rep(rownames(tab), times = ncol(tab)),
                         dir = rep(colnames(tab), each = nrow(tab)),
                         Number = as.vector(tab), stringsAsFactors = FALSE)
        p <- plotly::plot_ly(df, x = ~gene, y = ~Number, color = ~dir, type = "bar")
        plotly::layout(p, barmode = "stack", title = paste("Top", length(keep),
                                                           "genes enriched by DMPs"),
                       xaxis = list(title = ""), yaxis = list(title = "number of CpGs"),
                       margin = list(b = 140, t = 50))
    }

    volcano_plot <- function(sel) {
        if (!nrow(sel)) return(.champ_gui_empty("No significant DMP for these cutoffs."))
        d <- data.frame(logFC = sel$logFC,
                        neglog10P = -log10(pmax(sel$adj.P.Val, 1e-50)),
                        probe = rownames(sel), gene = sel$gene, stringsAsFactors = FALSE)
        n <- min(1500L, nrow(d))
        d <- d[order(d$neglog10P, decreasing = TRUE)[seq_len(n)], , drop = FALSE]
        p <- plotly::plot_ly(d, x = ~logFC, y = ~neglog10P, type = "scatter",
                             mode = "markers", text = ~paste0(probe, " (", gene, ")"),
                             marker = list(size = 5, opacity = 0.6, color = ~neglog10P,
                                           colorscale = "RdBu", reversescale = TRUE,
                                           showscale = FALSE))
        plotly::layout(p, title = paste0("Volcano plot for ", n, " significant CpGs"),
                       xaxis = list(title = "logFC"),
                       yaxis = list(title = "-log10 adjusted p"))
    }

    gene_plot <- function(sel) {
        if (!nrow(sel)) return(.champ_gui_empty("No significant CpG for this gene."))
        sel <- sel[order(sel$MAPINFO), , drop = FALSE]
        mat <- beta[rownames(sel), , drop = FALSE]
        xs <- seq_len(nrow(mat))
        pal <- .champ_gui_palette(n.gr)
        p <- plotly::plot_ly()
        for (i in seq_len(n.gr)) {
            g <- levels(pheno)[i]
            cols <- which(pheno == g)
            if (!length(cols)) next
            sub <- data.frame(
                x = rep(xs, times = length(cols)),
                Beta = as.vector(mat[, cols, drop = FALSE]),
                label = paste0(rep(rownames(mat), times = length(cols)), " | ",
                               rep(colnames(mat)[cols], each = nrow(mat))),
                stringsAsFactors = FALSE)
            p <- plotly::add_trace(p, data = sub, x = ~x, y = ~Beta, type = "scatter",
                                   mode = "markers", name = g, text = ~label,
                                   marker = list(size = 4, color = pal[i], opacity = 0.55))
            p <- plotly::add_trace(p, x = xs,
                                   y = rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE),
                                   type = "scatter", mode = "lines", name = paste(g, "mean"),
                                   line = list(color = pal[i], width = 3))
        }
        rng <- suppressWarnings(range(mat, na.rm = TRUE))
        if (!all(is.finite(rng))) rng <- c(0, 1)
        step <- max(0.02, diff(rng) * 0.08)
        p <- .champ_gui_add_track(p, c(
            .champ_gui_track(sel$feature, .champ_gui_feature_cols, rng[1] - 2 * step),
            .champ_gui_track(sel$cgi, .champ_gui_cgi_cols, rng[1] - 3.5 * step)))
        plotly::layout(p, title = paste0(sel$gene[1L], " (", nrow(sel), " CpG)"),
                       xaxis = list(title = "CpG in genomic order"),
                       yaxis = list(title = "beta", range = c(rng[1] - 4.5 * step, rng[2])),
                       margin = list(b = 60, t = 50))
    }

    cpg_plot <- function(cpgid) {
        cpgid <- trimws(as.character(cpgid)[1L])
        if (is.na(cpgid) || !nzchar(cpgid)) return(.champ_gui_empty("Enter a CpG ID."))
        if (!cpgid %in% rownames(beta)) {
            return(.champ_gui_empty(paste0("CpG '", cpgid, "' is not in 'beta'.")))
        }
        d <- data.frame(Beta = as.numeric(beta[cpgid, ]), Pheno = pheno,
                        Sample = colnames(beta), stringsAsFactors = FALSE)
        if (.champ_is_numeric(pheno)) {
            p <- plotly::plot_ly(d, x = ~Pheno, y = ~Beta, type = "scatter", mode = "markers",
                                 text = ~Sample, marker = list(size = 6, opacity = 0.6))
            p <- plotly::layout(p, xaxis = list(title = "phenotype"),
                                yaxis = list(title = "beta"))
        } else {
            p <- plotly::plot_ly(d, x = ~Pheno, y = ~Beta, type = "box",
                                 boxpoints = "all", jitter = 0.3, pointpos = 0,
                                 color = ~Pheno)
            p <- plotly::layout(p, xaxis = list(title = ""), yaxis = list(title = "beta"),
                                showlegend = FALSE)
        }
        plotly::layout(p, title = cpgid)
    }

    ui <- shiny::fluidPage(
        theme = .champ_gui_theme(),
        shiny::tags$head(shiny::tags$style(
            "#table, #genetable {font-size:75%;}",
            "#heatmap {height:70vh !important;}",
            "#geneplot {height:60vh !important;}",
            "#volcanoplot, #cpgplot, #cgibarplot, #featurebarplot,",
            "#featurecgibarplot, #geneenrich {height:45vh !important;}")),
        shiny::titlePanel("DMP Overview"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                width = 3,
                shiny::sliderInput("pvaluebin", "P value cutoff:", min = 0, max = 1,
                                   step = 0.025, value = 0.05),
                shiny::sliderInput("logFCbin", "abs(logFC) cutoff:", min = 0, max = 1,
                                   step = 0.025, value = 0),
                shiny::actionButton("go", "Submit"),
                shiny::hr(),
                shiny::textInput("genebin", "Gene symbol:", value = gene_default),
                shiny::actionButton("genesearch", "Submit"),
                shiny::hr(),
                shiny::textInput("cpgbin", "CpG ID:", value = rownames(DMP)[1L]),
                shiny::actionButton("cpgsearch", "Submit")),
            shiny::mainPanel(
                shiny::tabsetPanel(
                    shiny::tabPanel("DMP table",
                                    shiny::div(DT::dataTableOutput("table"),
                                               style = "font-size:75%")),
                    shiny::tabPanel("Heatmap", plotly::plotlyOutput("heatmap")),
                    shiny::tabPanel("Feature & CpG",
                                    shiny::fluidRow(
                                        shiny::column(6, plotly::plotlyOutput("cgibarplot")),
                                        shiny::column(6, plotly::plotlyOutput("featurebarplot")),
                                        shiny::column(12, plotly::plotlyOutput("featurecgibarplot")))),
                    shiny::tabPanel("Gene",
                                    plotly::plotlyOutput("geneplot"),
                                    shiny::div(DT::dataTableOutput("genetable"),
                                               style = "font-size:75%")),
                    shiny::tabPanel("Volcano",
                                    plotly::plotlyOutput("volcanoplot")),
                    shiny::tabPanel("CpG",
                                    plotly::plotlyOutput("geneenrich"),
                                    plotly::plotlyOutput("cpgplot"))
                )
            )
        )
    )

    server <- function(input, output, session) {
        cut_data <- shiny::eventReactive(input$go, {
            p <- as.numeric(input$pvaluebin)
            fc <- as.numeric(input$logFCbin)
            if (!is.finite(p)) p <- 1
            if (!is.finite(fc)) fc <- 0
            sel <- DMP[DMP$adj.P.Val <= p & abs(DMP$logFC) >= fc, , drop = FALSE]
            talk("  DMP.GUI: ", nrow(sel), " of ", nrow(DMP), " DMPs pass the cutoffs.")
            list(sel = sel,
                 cgi = prop_table(DMP$cgi, sel$cgi, sel$logFC, fc),
                 feature = prop_table(DMP$feature, sel$feature, sel$logFC, fc),
                 featcgi = prop_table(DMP$feat.cgi, sel$feat.cgi, sel$logFC, fc))
        }, ignoreNULL = FALSE)

        gene_data <- shiny::eventReactive(input$genesearch, {
            g <- trimws(input$genebin)
            p <- as.numeric(input$pvaluebin)
            fc <- as.numeric(input$logFCbin)
            if (!is.finite(p)) p <- 1
            if (!is.finite(fc)) fc <- 0
            if (!nzchar(g)) return(list(sel = DMP[0, , drop = FALSE], gene = g))
            sel <- DMP[grepl(g, DMP$gene, fixed = TRUE) & DMP$adj.P.Val <= p &
                       abs(DMP$logFC) >= fc, , drop = FALSE]
            list(sel = sel, gene = g)
        }, ignoreNULL = FALSE)

        cpg_data <- shiny::eventReactive(input$cpgsearch, {
            trimws(input$cpgbin)
        }, ignoreNULL = FALSE)

        output$table <- DT::renderDataTable(dmptable(cut_data()$sel))
        output$genetable <- DT::renderDataTable(dmptable(gene_data()$sel))
        output$heatmap <- plotly::renderPlotly(heat_plot(cut_data()$sel))
        output$cgibarplot <- plotly::renderPlotly(bar_plot(cut_data()$cgi, "CpG island"))
        output$featurebarplot <- plotly::renderPlotly(bar_plot(cut_data()$feature, "Feature"))
        output$featurecgibarplot <- plotly::renderPlotly(
            bar_plot(cut_data()$featcgi, "Feature - CpG island"))
        output$geneenrich <- plotly::renderPlotly(gene_enrich(cut_data()$sel))
        output$volcanoplot <- plotly::renderPlotly(volcano_plot(cut_data()$sel))
        output$geneplot <- plotly::renderPlotly(gene_plot(gene_data()$sel))
        output$cpgplot <- plotly::renderPlotly(cpg_plot(cpg_data()))
    }

    app <- shiny::shinyApp(ui = ui, server = server)
    talk("[<<<<<< DMP.GUI END >>>>>>>>>]")
    if (!isTRUE(runApp)) return(app)
    invisible(shiny::runApp(app))
}