# ChAMP2 — Chip Analysis Methylation Pipeline

[![R](https://img.shields.io/badge/R-4.2.3--4.6.2-blue)](https://cran.r-project.org)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.16--3.24-green)](https://bioconductor.org)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

ChAMP2 is a maintained, dependency-hardened fork of
[ChAMP](https://bioconductor.org/packages/ChAMP) (Tian et al., *Bioinformatics* 2017).
It performs quality control, normalisation, differential methylation, copy-number and
gene-set enrichment analysis for Illumina **450K**, **EPIC (v1)** and **EPICv2** arrays.

The upstream 2.29.1 source tree could only be installed on the R/Bioconductor release it
was written for. This fork installs and runs unchanged on **every R from 4.2.3 to 4.6.2**
(Bioconductor 3.16 – 3.24), and it tells you exactly what to do when an *optional*
feature is unavailable.

---

## 1. Quick start

```r
# 1. install every dependency (required + optional), any R 4.2.3 - 4.6.2
Rscript inst/scripts/install_deps.R
# ...or from R:
# source("https://raw.githubusercontent.com/hexiruofeng/ChAMP2/master/inst/scripts/install_deps.R")

# 2. install ChAMP itself
# R CMD INSTALL ChAMP2

library(ChAMP)
champ_check_env()                 # what is available, what is missing
```

One call runs the whole pipeline:

```r
myResults <- champ.process(directory = "IDAT/", arraytype = "EPICv2")
```

or step by step:

```r
myLoad  <- champ.load("IDAT/", arraytype = "450K")
myNorm  <- champ.norm(myLoad$beta, arraytype = "450K", method = "BMIQ", cores = 4)
myDMP   <- champ.DMP(myNorm, pheno = myLoad$pd$Sample_Group)
myDMR   <- champ.DMR(myNorm, pheno = myLoad$pd$Sample_Group, method = "DMRcate")
myGSEA  <- champ.GSEA(beta = myNorm, DMP = myDMP[[1]], DMR = myDMR, arraytype = "450K")
```

---

## 2. Compatibility

| R           | Bioconductor | Status     |
|-------------|--------------|------------|
| 4.2.3       | 3.16         | supported  |
| 4.3.x       | 3.17 / 3.18  | supported  |
| 4.4.x       | 3.19 / 3.20  | supported  |
| 4.5.x       | 3.21 / 3.22  | supported  |
| 4.6.0–4.6.2 | 3.23 / 3.24  | supported  |

Release 2.32.0 was validated end to end on **R 4.6.1 / Bioconductor 3.23**
(`R CMD check`: 0 errors, 0 warnings, 0 notes; 90 `testthat` assertions) and on
**R 4.3.3 / Bioconductor 3.18** (install + 90 `testthat` assertions). The
GitHub Actions matrix builds and checks R 4.2.3, 4.3.3, 4.4.3, 4.5.3 and 4.6.1.

The code avoids every construct whose behaviour changed in that window: no
`class(x) == "..."`, no `if (vector | vector)`, no `T`/`F`, no implicit
`stringsAsFactors`, no `aes_string()`, no `import()` of whole packages. Every access to an
optional package goes through `pkg::` and a guard that explains how to install it.

---

## 3. Dependencies: required, optional, and bundled

### Required (installed automatically by `install_deps.R`)

`ChAMPdata`, `limma`, `matrixStats`, `illuminaio`, `DNAcopy`,
`preprocessCore`, `bumphunter`, `sva`, `globaltest`,
`GenomicRanges`, `IRanges`, `S4Vectors`, `GenomeInfoDb`, `quadprog`.

### Optional (enable extra features)

| Package | Enables |
|---------|---------|
| `minfi` | `champ.load(method = "minfi")`, SWAN, FunctionalNormalization, control-probe SVD |
| `DMRcate` | the original DMRcate backend (a compatible implementation is bundled) |
| `goseq` + `geneLenDataBase` | `champ.GSEA(method = "goseq")` |
| `missMethyl` | `champ.GSEA(method = "gometh")` |
| `shiny`, `shinythemes`, `plotly`, `DT` | the five interactive GUIs |
| `dendextend` | coloured dendrograms in `champ.QC()` |
| `RPMM`, `impute`, `kpmt` | the historical backends (see below) |

### Bundled replacements

Several helper packages ChAMP used to import are no longer reliably installable. Their
algorithms are now part of ChAMP, are used automatically when the original package is
absent, and are tested against the originals when it is present:

| Removed/archived package | Replacement | Used by |
|---|---|---|
| `impute` (KNN) | `champ_knn_impute()` | `champ.load`, `champ.filter`, `champ.impute`, `champ.norm` |
| `RPMM` (`blc`) | `champ_blc()` | `champ.BMIQ()` |
| `kpmt` (archived on CRAN 2025-04-27) | `champ_kpmt()` | `champ.ebGSEA()` |
| `DMRcate` (not installable everywhere) | `champ_dmrcate()` | `champ.DMR(method = "DMRcate")` |
| `marray`, `wateRmelon`, `Illumina450ProbeVariants.db` | no longer needed (the imports were unused) | — |

When the original package *is* installed, it is used by default, so results are identical
to the legacy pipeline. Set `use_impute = FALSE` (or `dmrcate.backend = "internal"`)
to force the bundled implementation.

---

## 4. What changed relative to 2.29.1

* **Installation**: `Depends: R (>= 3.3), minfi, DMRcate, DT, RPMM` is gone; only
  `R (>= 4.2.3)` plus 14 required packages. `minfi` is now optional.
* **Annotations**: all `data()` calls were replaced by one cached loader, because several
  objects changed file names between ChAMPdata releases (`probe.features.epicv1` never
  existed as a file, `bloodCtl` lives in `champBloodCtls.rda`, `illumina450Gr` has been
  removed). EPICv1/EPICv2 now work end to end.
* **ProbeLasso** builds its GRanges from the annotation instead of requiring the removed
  `illumina450Gr` object.
* **champ.process** no longer calls the non-existent `champ.EpiMod()`, and a failing step
  no longer aborts the whole pipeline.
* **Bugs fixed**: the `stringsAsFactor` typo in `champ.import`, the undefined `h2` in
  `DMP.GUI`, `sample(1:n, nfit)` crashing for small arrays in BMIQ, `if (x | y)`
  conditions, `else (stop(...))` blocks, undefined `RemainSample` in `champ.filter`,
  `champ.SVD` failing on integer matrices, `champ.CNA` failing without control probes.
* **Parallelism** uses the base `parallel` package, so `cores` works on Windows too.

See [NEWS.md](NEWS.md) for the full list, [COMPATIBILITY.md](COMPATIBILITY.md) for the
design notes, and [inst/INSTALL-zh.md](inst/INSTALL-zh.md) for a Chinese
installation guide.

---

## 5. The GUIs

```r
DMP.GUI(myDMP[[1]], beta = myNorm, pheno = myLoad$pd$Sample_Group)
DMR.GUI(myDMR[[1]], beta = myNorm, pheno = myLoad$pd$Sample_Group)
Block.GUI(myBlock$Block, beta = myNorm, pheno = myLoad$pd$Sample_Group)
QC.GUI(beta = myNorm, pheno = myLoad$pd$Sample_Group)
CpG.GUI(CpG = rownames(myDMP[[1]])[1:100])
```

They require `shiny`, `plotly`, `DT` and `shinythemes`; pass `runApp = FALSE`
to get the app object without launching it (useful in tests).

---

## 6. Troubleshooting

| Symptom | Fix |
|---|---|
| `there is no package called 'minfi'` | `ChAMP::champ_install_deps()`, or use `champ.load(method = "ChAMP")` which does not need minfi |
| `No probe annotation available` | reinstall ChAMPdata: `BiocManager::install("ChAMPdata", force = TRUE)` |
| `Could not load dataset 'bloodCtl'` | your ChAMPdata is too old; `ChAMP::champ_install_deps()` |
| BMIQ fails for one sample | low-quality sample: check `champ.QC()`, the sample may lack a beta distribution |
| `champ.DMP()` finds nothing | relax `adjPVal` or check the group labels |
| Windows: `cores > 1` is slow | BMIQ parallelises with PSOCK clusters; use `cores = 1` on small data |

---

## 7. 中文说明

ChAMP2 是对 Bioconductor 上 [ChAMP](https://bioconductor.org/packages/ChAMP) 的加固分支，
目标只有一个：**同一份源码，在 R 4.2.3 到 R 4.6.2（Bioconductor 3.16–3.24）上都能安装、都能跑通**。

* **一键装依赖**：`Rscript inst/scripts/install_deps.R`（或 `ChAMP::champ_install_deps()`），
  自动识别 R/Bioconductor 版本、优先使用二进制包，必装与选装分开，装不上的选装包不会阻塞主流程。
* **不再强制依赖 minfi / DMRcate / RPMM / impute / marray**：这些包在部分 R 版本或系统上装不上，
  现在改为可选；对应算法（KNN 填补、BMIQ 的 beta 混合模型、KPMT、DMRcate 区域识别）已内置，
  装上原包时仍优先调用原包，结果与老版本一致。
* **注释数据**：ChAMPdata 里若干数据集改过名（`probe.features.epicv1` 实际不存在、
  `bloodCtl` 在 `champBloodCtls.rda` 里、`illumina450Gr` 已删除），本版本用统一的
  缓存加载器处理，EPICv1 / EPICv2 全流程可用。
* **修掉的老 bug**：`stringsAsFactor` 拼写、`champ.process` 调用不存在的 `champ.EpiMod`、
  BMIQ 在小阵列上 `sample()` 越界、`if (x | y)` 向量条件、`champ.filter` 里未定义的
  `RemainSample`、`DMP.GUI` 里未定义的 `h2` 等。
* **Windows 也能多核**：并行改用 base `parallel`，不再依赖 `doParallel` 的 fork 语义。

安装：

```bash
git clone https://github.com/hexiruofeng/ChAMP2.git
Rscript ChAMP2/inst/scripts/install_deps.R
R CMD INSTALL ChAMP2
```

```r
library(ChAMP)
champ_check_env()                                   # 查看依赖状态
res <- champ.process("IDAT目录", arraytype = "EPICv2")  # 一条命令跑完整流程
```

---

## 8. Citation

Yuan Tian, Tiffany J Morris, Amy P Webster, Zhen Yang, Stephan Beck, Andrew Feber,
Andrew E Teschendorff. *ChAMP: updated methylation analysis pipeline for Illumina
BeadChips.* Bioinformatics (2017) btx513.

The bundled algorithms:
* BMIQ — Teschendorff et al., *Bioinformatics* (2013).
* DMRcate — Peters et al., *Epigenetics & Chromatin* (2015).
* ComBat — Johnson, Li & Rabinovic, *Biostatistics* (2007).
* KPMT — Conover, *Practical Nonparametric Statistics*.

## 9. License

GPL-3, as the original ChAMP package.