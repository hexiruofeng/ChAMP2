# ChAMP2 2.32.0

## Compatibility (the reason for this release)

* Supports **R 4.2.3 – 4.6.2** and **Bioconductor 3.16 – 3.24** from a single
  source tree.
* DESCRIPTION now requires only `R (>= 4.2.3)`. The former hard dependencies on
  `minfi`, `DMRcate`, `DT`, `RPMM`, `impute`, `marray`, `wateRmelon` and
  `Illumina450ProbeVariants.db` are gone; `minfi` and the GUI packages moved to
  `Suggests` and every use of them is guarded with an actionable message.
* Removed all `import()` directives in favour of `pkg::` calls plus a small
  `importFrom` block, so namespace clashes (for example `plotly::last_plot` vs
  `ggplot2::last_plot`) can no longer break the package at load time.

## Bundled replacements for unreachable dependencies

* `champ_knn_impute()` — KNN imputation, drop-in for `impute::impute.knn()`
  (same arguments: `k`, `rowmax`, `colmax`, `maxp`). Used automatically when
  `impute` is not installed; `impute` is still preferred when present.
* `champ_blc()` — beta-mixture EM used by BMIQ, re-implementation of
  `RPMM::blc()` (archived on CRAN). Same return components (`a`, `b`, `eta`,
  `mu`, `w`, `llike`).
* `champ_kpmt()` — exact known-population median test, replacement for
  `kpmt::kpmt()` (archived on CRAN 2025-04-27). Returns the same seven-column
  data frame, so `[[4]]` (the p-value) still works.
* `champ_dmrcate()` — DMRcate-compatible DMR calling (kernel chi-square test,
  Stouffer/HM/Fisher combination, identical output columns). Used when
  `DMRcate` is not installed or cannot handle the array (`dmrcate.backend`).

## Fixed

* `champ.import()`: `stringsAsFactor` typo (`read.csv` silently used the
  default); detection p-values are now computed against the negative controls
  even when their names changed; zero-variance control probes no longer produce
  `NaN`; duplicate CpG identifiers are collapsed.
* `champ.filter()`: `RemainSample` was undefined when `filterDetP = FALSE`;
  `else (stop(...))` constructs replaced by proper branches; filtering now
  validates dimensions instead of silently recycling.
* `champ.process()`: called the removed `champ.EpiMod()` (now
  `champ.ebGSEA()`); a failing step aborted the whole pipeline (now reported and
  skipped unless `strict = TRUE`); the CNA step ignored `resultsDir`.
* `champ.BMIQ()`: `sample(1:length(beta), nfit)` failed whenever a sample had
  fewer than `nfit` probes; plots are written to `plotDir` instead of changing
  the working directory; constant type-I/II subsets now give a clear error.
* `champ.norm()`: `%dopar%` exported `blc` even when RPMM was absent; the
  working directory was changed globally; NA values are now imputed before
  normalisation.
* `champ.SVD()`: failed for integer matrices (`is.data.frame(myNorm)` tested the
  wrong variable), the heatmap legend had swapped coordinates, and `pd` columns
  with a single value crashed the loop.
* `champ.DMP()` / `champ.DMR()` / `champ.GSEA()`: `class(x) == "numeric"` and
  `if (a | b)` replaced; multi-group comparisons no longer drop the group names;
  the numeric-covariate branch reports the correct coefficient.
* `champ.DMR(method = "ProbeLasso")`: used the `illumina450Gr` /
  `illuminaEPICGr` data objects that ChAMPdata no longer ships; the GRanges are
  now built from the probe annotation, so ProbeLasso works for 450K, EPICv1 and
  EPICv2.
* `champ.CNA()`: crashed when the array lacked negative controls; chromosome
  levels are normalised; NA intensities no longer abort the segmentation.
* `champ.refbase()`: `cellFrac[,-1*which.min(...)]` produced a vector for small
  reference sets; the design matrix is now built explicitly.
* `champ.runCombat()`: `apply(pd, 2, ...)` turned every column into character;
  a batch variable identical to the phenotype is rejected with a clear message.
* GUIs: undefined `h2` in `DMP.GUI`, `ggplot2::diamonds`/`economics` used as if
  they were base datasets, `plotly` NSE calls that broke on plotly >= 4.10.

## Changed

* `arraytype` accepts every spelling ("450K", "EPIC", "EPICv1", "EPICv2",
  "EPIC-1", "HM450", ...); `"EPIC"` still means EPICv2 (as since 2.29) but now
  warns about it once.
* Annotation access goes through a cached loader that understands object names
  that differ from file names (`bloodCtl` in `champBloodCtls.rda`), renamed
  datasets (`probe.features.epicv1`) and missing datasets (`illumina450Gr`),
  with minfi annotation as a fallback.
* `arraytype` is normalised in every function; chromosome and CpG-island labels
  are canonicalised ("chr1" -> "1", "OpenSea" -> "opensea") so that the same
  code works with ChAMPdata 2.6 and 2.31.
* `champ.QC()` no longer needs minfi: MDS, density and dendrogram plots are
  implemented locally (dendextend is used when available for coloured leaves).
* Parallel execution uses base `parallel` (PSOCK), so `cores > 1` works on
  Windows as well.
* `champ.install_deps()` and `inst/scripts/install_deps.R` install the whole
  dependency set, required and optional, for the running R/Bioconductor version.
* Every function gained a `verbose` argument; set `options(ChAMP.quiet = TRUE)`
  to silence everything.

## New exported functions

`champ_knn_impute()`, `champ_blc()`, `champ_kpmt()`, `champ_dmrcate()`,
`champ_check_env()`, `champ_install_deps()`, `champ_dependencies()`,
`champ.read()`, `champ.colname()`.

## Tests

A `testthat` suite covering the compatibility layer (KNN imputation, beta
mixture EM, exact KPMT against exhaustive enumeration, DMRcate region calling),
the filtering engine and the annotation loader.
