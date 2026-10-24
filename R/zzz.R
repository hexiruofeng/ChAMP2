# =============================================================================
#  Package startup
# =============================================================================

utils::globalVariables(c(
    "pve", "names", "probe", "proportion", "cgi", "gene", "Number", "Variance",
    "Value", "ID", "Sample", "logFC", "adj.P.Val", "P.Value", "Beta_Value",
    "Pheno", "size", "info", "x", "y", "pos", "Mean", "Fitted", "cut",
    "showtext", "dmr", "chrom", "start", "end", "nREP", "AUC", "P", "adjP",
    "chr", "MAPINFO", "CHR", "mean", "sd", "group", "pheno", "beta", "M",
    "probeID", "TargetID", "MASK_general", "Type", "feature", "feat.cgi",
    "seqnames", "width", "strand", "no.cpgs", "min_smoothed_fdr", "Stouffer",
    "HMFDR", "Fisher", "maxdiff", "meandiff", "overlapping.genes",
    "sampleNames<-", "mcols<-", "seqlevels<-", "seqlevels", "values<-"
))

.onAttach <- function(libname, pkgname) {
    mes <- paste0(
        "       ___ _      _   __  __ ___ \n",
        "      / __| |_   /_\\ |  \\/  | _ \\\n",
        "     | (__| ' \\ / _ \\| |\\/| |  _/\n",
        "      \\___|_||_/_/ \\_\\_|  |_|_|  \n",
        "      ------------------------------\n")
    ver <- tryCatch(as.character(utils::packageVersion("ChAMP")),
                    error = function(e) "2.32.0")
    packageStartupMessage(
        ">> ChAMP (ChAMP2) ", ver, " loaded <<\n", mes,
        "    R ", getRversion(), " | Bioconductor ",
        tryCatch(as.character(BiocManager::version()), error = function(e) "unknown"),
        "\n    If you have any question or suggestion about ChAMP, please email champ450k@gmail.com.\n",
        "    Please cite:\n",
        "    Yuan Tian, Tiffany J Morris, Amy P Webster, Zhen Yang, Stephan Beck, Andrew Feber,\n",
        "    Andrew E Teschendorff; ChAMP: updated methylation analysis pipeline for Illumina\n",
        "    BeadChips, Bioinformatics, btx513, https://doi.org/10.1093/bioinformatics/btx513\n",
        "    Run ChAMP::champ_check_env() to inspect the optional dependencies.")
}

.onLoad <- function(libname, pkgname) {
    ## Respect the user's locale without changing it globally.
    invisible(NULL)
}
