# Compatibility notes

This document explains *why* ChAMP 2.32 looks the way it does, so that the next
person who has to touch it does not re-introduce the problems.

## The support window

| R | Bioconductor | Notes |
|---|--------------|-------|
| 4.2.3 | 3.16 | oldest supported; base `|>`, `gregexpr()` with `useBytes` |
| 4.3.x | 3.17, 3.18 | S4Vectors 0.38 changes `DataFrame` subsetting |
| 4.4.x | 3.19, 3.20 | `minfi` moved to a new annotation layout |
| 4.5.x | 3.21, 3.22 | matrixStats 1.0 deprecates `rowMads`/`colMads` |
| 4.6.0–4.6.2 | 3.23, 3.24 | current release; `DNAratio`/EPICv2 default |

A single source tree has to compile and pass checks on all of them, which rules
out anything that is only available in one release.

## Rules the code follows

1. **No whole-package imports.** `import(pkg)` exports every symbol of a package
   into the ChAMP namespace. Between Bioconductor releases these sets change
   (for example `plotly::last_plot` disappeared, `matrixStats` renamed
   functions), and a package that imports everything breaks at load time with
   "object not found". Every call is `pkg::fun()`; the NAMESPACE only imports
   the handful of functions that are used unqualified.
2. **Optional packages are optional.** A feature that needs `minfi` checks for it
   with `.champ_require("minfi", "champ.load(method = 'minfi')")` and stops with
   a message naming the function, the package and the install command. Nothing
   in `.onLoad` touches an optional package.
3. **No behaviour that changed in the window.** `class(x) == "numeric"`,
   `if (length(x) > 1)`, `if (a | b)`, `T`/`F`, `stringsAsFactors`,
   `data.frame(..., stringsAsFactors = TRUE)`, `aes_string()` and
   `matrixStats::rowMads()` are all gone.
4. **Annotation is data, not code.** The historical code called `data(probe.features)`
   in twenty places with different expectations. Now `.champ_probe_features(arraytype)`
   is the only entry point: it tries the ChAMPdata object names in order, falls
   back to a direct `load()` of the `.rda` file, then to the minfi annotation,
   canonicalises the columns (`CHR` `"chr1"` -> `"1"`, `cgi` `"OpenSea"` ->
   `"opensea"`, `Type` `"II"`) and caches the result.
5. **Deterministic randomness.** Set seeds are local (`set.seed()` inside the
   function that needs it) and `champ_knn_impute()` seeds its row sample.

## The bundled algorithms

### `champ_knn_impute()` (replaces `impute::impute.knn()`)

The original C implementation computes, for every missing entry of a probe, the
`k` nearest *samples* (columns) among those with an observed value, and averages
them. The R implementation reproduces that:

* probes with more than `rowmax` missing values are filled with column means
  (this is what `impute.knn()` does, with a warning);
* samples with more than `colmax` missing values are rejected, or mean-filled
  when `strict = FALSE` (ChAMP uses `strict = FALSE` so that a bad sample cannot
  abort a whole run);
* sample-to-sample distances come from at most `maxp` complete probes
  (`impute.knn()` splits the matrix into two clusters instead; for ChAMP's
  purpose - filling a handful of failed probes - the results agree to within
  the noise of the imputation itself);
* rows sharing a missingness pattern share the neighbour search, which keeps the
  cost at `O(#patterns x #missing columns x k)` instead of `O(#missing values)`.

### `champ_blc()` (replaces `RPMM::blc()`)

RPMM was archived on CRAN. BMIQ only uses the intercept-only, single-covariate
case, so the class model reduces to a finite beta mixture:

* E-step: `w_ik proportional to eta_k * dbeta(y_i; a_k, b_k)` with a
  log-sum-exp stabilisation;
* M-step: `a_k, b_k` by weighted maximum likelihood (`stats::optim`, BFGS,
  started from the method-of-moments estimate), `eta_k = mean(w_k)`;
* convergence on `max|mu_k - mu_k^old| < tol`, exactly like `blc()`;
* the `mu`, `a`, `b`, `eta`, `w`, `llike` components have the same shapes, so
  `champ.BMIQ()` needed no change beyond the call.

### `champ_kpmt()` (replaces `kpmt::kpmt()`)

The known-population median test is an exact combinatorial test. For a sample of
size `s` from `N` ordered values with minimal median `w`:

  P(min-median <= w) = s!(N-s)!/(N!(m-1)!(s-m)!) * sum over x in U of
                       (x-1)!(N-x)! / ((x-m)!(N-x-s+m)!)

with `m = ceiling(s/2)` and `U` the ranks that are below `w` and leave enough
room for the sample. `lgamma()` makes the computation stable and the sum is
vectorised, so no package is needed. The returned data frame has the same seven
columns as `kpmt::kpmt()`, including `p.value` in position four.

### `champ_dmrcate()` (DMRcate-compatible)

The DMRcate algorithm (Peters et al. 2015) as implemented here:

1. moderated t-statistics per probe (limma), `X2 = t^2`;
2. per chromosome, local kernel sums of `X2` with a Gaussian kernel of bandwidth
   `sigma = lambda / C` (default 500 bp), then Satterthwaite's approximation:
   `a = skk/sk`, `b = sk^2/skk`, smoothed `X2 = sky/a`, `p = pchisq(., b,
   lower.tail = FALSE)`;
3. probes are "significant" when their limma FDR is below `fdr`; consecutive
   significant probes (gap <= `lambda`) form a region;
4. region statistics: `min_smoothed_fdr`, Stouffer (`qnorm` mean), harmonic mean,
   Fisher combination, `maxdiff`, `meandiff`, adjusted with BH;
5. regions shorter than `min.cpgs` are dropped and the table is sorted by
   `min_smoothed_fdr` and `no.cpgs`.

The result has the same columns as `as.data.frame(DMRcate::extractRanges())`
except `overlapping.genes`, which is computed from the ChAMP annotation instead
of ExperimentHub (so no network access is required).

## What is *not* supported

* `minfi`'s `preprocessFunnorm()` on EPICv2 requires the matching annotation
  package; ChAMP guards it and reports the missing package.
* The historical `champ.EpiMod()` (EpiMod/HPRD) was removed upstream long ago;
  `champ.process(runEpiMod = TRUE)` now runs `champ.ebGSEA()` under that name.
* 27K arrays are accepted by `champ.refbase()` only; the rest of the pipeline
  targets 450K/EPIC/EPICv2.
