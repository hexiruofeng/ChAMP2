# ChAMP2 安装与使用说明（中文）

本说明对应 **ChAMP 2.32.0（ChAMP2 分支）**，支持 **R 4.2.3 – 4.6.2**
（Bioconductor 3.16 – 3.24）。

---

## 一、准备 R

先安装 R（建议 4.6.x，但 4.2.3 以上都能用）：

* Windows / macOS：<https://cran.r-project.org/>
* Linux：用系统包管理器安装 `r-base`、`r-base-dev`、`libcurl4-openssl-dev`、
  `libxml2-dev`、`libssl-dev`、`libpng-dev` 等编译依赖。

装好后在命令行确认：

```bash
R --version
Rscript --version
```

> Windows 用户如果之前装过 R，请确认 `Rscript.exe` 在 PATH 中；
> 也可以在 R 里直接执行后面所有步骤，无需命令行。

---

## 二、安装依赖（一步完成）

ChAMP 需要 14 个必装包和若干可选包。仓库自带一个安装脚本，会自动识别
R / Bioconductor 版本并优先使用二进制包：

```bash
Rscript ChAMP2/inst/scripts/install_deps.R
```

常用参数：

| 参数 | 作用 |
|---|---|
| `--required` | 只装必装包（最快） |
| `--check` | 只检查，不安装；缺必装包时退出码为 1 |
| `--source` | 强制源码编译（需要 Rtools / 编译工具链） |
| `--lib=D:/Rlibs` | 装到指定库目录 |
| `--quiet` | 少打印 |

R 里也可以直接调用：

```r
install.packages("ChAMP_2.32.0.tar.gz", repos = NULL, type = "source")  # 先装 ChAMP
ChAMP::champ_install_deps()      # 再补齐依赖
ChAMP::champ_check_env()         # 查看状态
```

如果网络访问 CRAN/Bioconductor 慢，可以先把 CRAN 镜像换成国内源：

```r
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
```

---

## 三、安装 ChAMP 包

```bash
# 方式一：命令行
R CMD INSTALL ChAMP_2.32.0.tar.gz

# 方式二：R 内
install.packages("ChAMP_2.32.0.tar.gz", repos = NULL, type = "source")

# 方式三：源码目录
R CMD INSTALL ChAMP2
```

安装完成后：

```r
library(ChAMP)
champ_check_env()
```

输出中 `ok` 表示已安装；`MISSING` 只可能出现在必装包上，按提示运行
`ChAMP::champ_install_deps()` 即可；`optional` 表示选装，不影响主流程。

---

## 四、快速跑通一个流程

```r
library(ChAMP)

# 一条命令跑完（载入→过滤→QC→归一化→SVD→ComBat→DMP→DMR→Block→GSEA→CNA）
res <- champ.process(directory = "IDAT目录", arraytype = "EPICv2")

# 或者分步
myLoad <- champ.load("IDAT目录", arraytype = "450K")
myNorm <- champ.norm(myLoad$beta, method = "BMIQ", arraytype = "450K", cores = 4)
myDMP  <- champ.DMP(myNorm, pheno = myLoad$pd$Sample_Group)
myDMR  <- champ.DMR(myNorm, pheno = myLoad$pd$Sample_Group, method = "DMRcate")
myGSEA <- champ.GSEA(beta = myNorm, DMP = myDMP[[1]], DMR = myDMR)
```

图形界面：

```r
DMP.GUI(myDMP[[1]], beta = myNorm, pheno = myLoad$pd$Sample_Group)
QC.GUI(beta = myNorm, pheno = myLoad$pd$Sample_Group)
CpG.GUI(rownames(myDMP[[1]])[1:100])
```

---

## 五、和旧版 ChAMP 2.29.1 的差别（重要）

1. **依赖变轻**：`minfi`、`DMRcate`、`DT`、`RPMM`、`impute`、`marray`、
   `wateRmelon` 都不再是硬依赖。装不上不会导致 `library(ChAMP)` 失败。
2. **内置算法**：
   * `champ_knn_impute()` —— 替代 `impute::impute.knn()`（KNN 填补）
   * `champ_blc()` —— 替代 `RPMM::blc()`（BMIQ 的 beta 混合模型）
   * `champ_kpmt()` —— 替代已下架的 `kpmt`（精确的已知总体中位数检验）
   * `champ_dmrcate()` —— DMRcate 算法的等价实现（无需 ExperimentHub 联网）
   
   装了原包时默认仍调用原包，结果与旧版一致；想强制用内置实现，可设
   `use_impute = FALSE` 或 `dmrcate.backend = "internal"`。
3. **注释数据统一读取**：ChAMPdata 里 `probe.features.epicv1` 实际不存在、
   `bloodCtl` 在 `champBloodCtls.rda` 中、`illumina450Gr` 已被删除——新代码
   会自动按多个候选名查找并做列名规范化，EPICv1 / EPICv2 全流程可用。
4. **修掉的老问题**：`stringsAsFactor` 拼写错误、`champ.process()` 调用不存在的
   `champ.EpiMod()`、BMIQ 在小阵列上 `sample()` 越界、`if (x | y)` 向量条件、
   `champ.filter()` 中未定义的 `RemainSample`、`DMP.GUI()` 中未定义的 `h2`、
   `champ.CNA()` 无对照探针时崩溃等。
5. **并行**：改用 base `parallel`（PSOCK），Windows 上 `cores > 1` 也能用。

---

## 六、常见问题

**Q1：`there is no package called 'minfi'`**
A：`ChAMP::champ_install_deps()`；或改用 `champ.load(method = "ChAMP")`，
   这条路不依赖 minfi（但 SWAN / FunctionalNormalization 仍需要它）。

**Q2：`Could not load dataset 'bloodCtl'` / 注释为空**
A：ChAMPdata 太旧或损坏。重装：`BiocManager::install("ChAMPdata", force = TRUE)`。

**Q3：`kpmt` 装不上（CRAN 已下架）**
A：正常。ChAMP 内置 `champ_kpmt()`，`champ.ebGSEA()` 的 `P(KPMT)` 列照常输出。

**Q4：R 4.6 下 `tabulate(..., weights=)` 报错**
A：R 4.6 移除了这个参数，本版本已改为内部的 `.champ_weighted_sum()`，不会触发。

**Q5：BMIQ 某个样本报错**
A：该样本 beta 分布异常（质控差）。先用 `champ.QC()` 看密度图，必要时剔除；
   也可改用 `method = "PBC"`。

**Q6：想完全离线安装**
A：在有网机器上运行 `Rscript tools/bundle_deps.R --dest=D:/champ2-deps` 下载全部
   依赖二进制包，拷贝到离线机器后运行
   `Rscript tools/bundle_deps.R --install --dest=D:/champ2-deps`。

---

## 七、自检

```bash
# 不依赖任何第三方包，检查源码与内置算法（适合装包前先验证）
Rscript tools/smoke_test.R

# 完整单元测试
Rscript tests/testthat.R
```

本项目在 R 4.6.1 / Bioconductor 3.23 上通过 `R CMD check`（0 error / 0 warning /
0 note）与 79 项单元测试，在 R 4.3.3 上通过源码解析与全部内置算法自检。
