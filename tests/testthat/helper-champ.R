# Shared helpers for the ChAMP test suite -------------------------------

skip_if_no_champdata <- function() {
    testthat::skip_if_not_installed("ChAMPdata")
}

# A small, deterministic beta matrix using real 450K probe identifiers so the
# annotation lookups behave exactly as in a real analysis.
make_beta <- function(n_probe = 500L, n_sample = 8L, seed = 42L) {
    skip_if_no_champdata()
    pf <- ChAMP:::.champ_probe_features("450K")
    ids <- rownames(pf)[!is.na(pf$MAPINFO)][seq_len(n_probe)]
    set.seed(seed)
    b <- matrix(stats::rbeta(length(ids) * n_sample, 2, 5),
                nrow = length(ids), ncol = n_sample,
                dimnames = list(ids, paste0("S", seq_len(n_sample))))
    ## add a group effect on the first 50 probes
    grp <- rep(c("A", "B"), each = n_sample / 2)
    b[seq_len(50), grp == "B"] <- pmin(1 - 1e-6, b[seq_len(50), grp == "B"] + 0.25)
    attr(b, "group") <- grp
    b
}

quiet_champ <- function(expr) {
    old <- options(ChAMP.quiet = TRUE)
    on.exit(options(old), add = TRUE)
    suppressMessages(expr)
}

clear_champ_warnings <- function() {
    env <- ChAMP:::.champ_state$warned
    if (length(ls(env))) rm(list = ls(env), envir = env)
    invisible(NULL)
}
