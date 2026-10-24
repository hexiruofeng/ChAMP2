# Dependency snapshot

`dependency-versions-R4.6-Bioc3.23.csv` records the exact package versions used
for the validation run of ChAMP 2.32.0:

* R 4.6.1 (2026-06-24 ucrt)
* Bioconductor 3.23
* ChAMPdata 2.31.1 (GitHub/YuanTian1991, the version ChAMP 2.29 was written for;
  Bioconductor 2.44 also works with the same loader)
* listed by `ChAMP::champ_check_env()`

Other R/Bioconductor combinations pick their own versions automatically:

| R | Bioconductor | install with |
|---|---|---|
| 4.2.3 | 3.16 | `Rscript inst/scripts/install_deps.R` |
| 4.3.x | 3.17 / 3.18 | `Rscript inst/scripts/install_deps.R` |
| 4.4.x | 3.19 / 3.20 | `Rscript inst/scripts/install_deps.R` |
| 4.5.x | 3.21 / 3.22 | `Rscript inst/scripts/install_deps.R` |
| 4.6.x | 3.23 / 3.24 | `Rscript inst/scripts/install_deps.R` |

To reproduce a snapshot on another machine, install exactly these versions with
`remotes::install_version()` / `BiocManager::install(version = "...")`, or run
`ChAMP::champ_check_env()` to see what is missing.

For an offline installation run `Rscript tools/bundle_deps.R --dest=<dir>` on a
connected machine and `--install` on the target machine.
