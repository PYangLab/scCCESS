# scCCESS <a href="https://github.com/PYangLab/scCCESS"><img src="https://i.imgur.com/frCIqJk.png" title="SnapCCESS hex sticker" align="right" height="138" /></a>

scCCESS implements **single-cell Consensus Clusters of Encoded Subspaces** for
single-cell RNA-seq clustering and cell-type number estimation. It combines
random gene subspaces, autoencoder-based low-dimensional encodings, and
consensus clustering to produce robust cell group assignments.

The package accompanies two papers:

- Yu, L., Cao, Y., Yang, J. Y. H. & Yang, P. Benchmarking clustering algorithms
  on estimating the number of cell types from single-cell RNA-sequencing data.
  *Genome Biology* 23, 49 (2022). <https://doi.org/10.1186/s13059-022-02622-0>
- Geddes, T. A., Kim, T., Nan, L., Burchfield, J. G., Yang, J. Y. H., Tao, D. &
  Yang, P. Autoencoder-based cluster ensembles for single-cell RNA-seq data
  analysis. *BMC Bioinformatics* 20, 660 (2019).
  <https://doi.org/10.1186/s12859-019-3179-5>

## Installation

Install the development version from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("PYangLab/scCCESS")
```

scCCESS uses Keras and TensorFlow for autoencoder training. For a CPU-only
setup, install the R packages and TensorFlow backend with:

```r
install.packages(c("keras", "tensorflow", "clue", "aricode"))
tensorflow::install_tensorflow()
```

GPU acceleration is optional and depends on a local TensorFlow/CUDA or
TensorFlow/ROCm installation. If TensorFlow is already managed through Conda,
Python, or a cluster module, configure the R `tensorflow` package to use that
environment before running scCCESS.

## Quick Start

```r
library(scCCESS)
library(SingleCellExperiment)

data("sce", package = "scCCESS")
dat <- SingleCellExperiment::counts(sce)

dat <- prefilter(dat)

k_result <- estimate_k(
  dat,
  seed = 1,
  cluster_func = function(x, centers) {
    set.seed(42)
    kmeans(x, centers)
  },
  criteria_method = "NMI",
  krange = 5:15,
  ensemble_sizes = 10,
  cores = 1
)

clusters <- ensemble_cluster(
  dat,
  seed = 1,
  cluster_func = function(x) {
    set.seed(1)
    kmeans(x, centers = k_result$ngroups)
  },
  genes_as_rows = TRUE,
  ensemble_sizes = 10,
  cores = 1,
  verbose = 0,
  scale = FALSE,
  batch_size = 64
)
```

## Input Data

Most functions accept a matrix-like object containing scRNA-seq expression
values. By default, `encode()` and `ensemble_cluster()` assume cells are rows and
genes are columns. Set `genes_as_rows = TRUE` when genes are rows and cells are
columns, which is the orientation returned by many Bioconductor count matrices.

`prefilter()` expects raw counts with genes as rows and cells as columns. It
filters low-quality cells, optionally removes genes with zero counts, and returns
log2-transformed counts.

Missing values are not supported. Replace or remove `NA` values before running
the package.

## Main Functions

| Function | Purpose |
| --- | --- |
| `prefilter()` | Basic count and gene filtering before clustering. |
| `encode()` | Train an autoencoder on a random gene subspace and return cell embeddings. |
| `ensemble_cluster()` | Build consensus clusters from repeated encoded subspaces. |
| `estimate_k()` | Estimate the number of cell types with clustering stability metrics. |

## Using SIMLR

scCCESS can use any clustering function that accepts a cell-by-feature matrix.
For SIMLR, `SIMLR_Large_Scale()` is usually preferred for speed:

```r
library(SIMLR)

k_result <- estimate_k(
  dat,
  seed = 1,
  cluster_func = function(x, centers) {
    set.seed(42)
    SIMLR_Large_Scale(t(x), c = centers, kk = 15)
  },
  criteria_method = "NMI",
  krange = 5:15,
  ensemble_sizes = 10,
  cores = 1
)
```

Some users have reported installation/runtime issues with the official SIMLR
package. If needed, an alternative fork without the optional t-SNE step has been
used previously:

```r
devtools::install_github("yulijia/SIMLR", ref = "master")
```

## Citation

If you use scCCESS, please cite the relevant method papers:

```bibtex
@article{Yu2022scCCESS,
  title = {Benchmarking clustering algorithms on estimating the number of cell types from single-cell RNA-sequencing data},
  author = {Yu, Lijia and Cao, Yue and Yang, Jean Yee Hwa and Yang, Pengyi},
  journal = {Genome Biology},
  year = {2022},
  volume = {23},
  number = {49},
  doi = {10.1186/s13059-022-02622-0},
  url = {https://doi.org/10.1186/s13059-022-02622-0}
}

@article{Geddes2019scCCESS,
  title = {Autoencoder-based cluster ensembles for single-cell RNA-seq data analysis},
  author = {Geddes, Thomas A. and Kim, Taiyun and Nan, Lihao and Burchfield, James G. and Yang, Jean Yee Hwa and Tao, Dacheng and Yang, Pengyi},
  journal = {BMC Bioinformatics},
  year = {2019},
  volume = {20},
  number = {660},
  doi = {10.1186/s12859-019-3179-5},
  url = {https://doi.org/10.1186/s12859-019-3179-5}
}
```

## License

scCCESS is distributed under the GPL-3 license.
