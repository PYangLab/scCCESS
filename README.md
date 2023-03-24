# scCCESS for cell clustering and number of cell type estimation from scRNA-seq data

scCCESS stands for single-cell Consensus Clusters of Encoded Subspaces. This repository stores the R implementation of scCCESS with stability-based approach for estimating number of cell types.


## Installation

```
if (!require("devtools")) install.packages("devtools")
devtools::install_github('PYangLab/scCCESS')
```

Package required to use this package:

- keras (available on CRAN) See https://keras.rstudio.com/reference/install_keras.html for more details - this package can be sped up substantially if the user has a supported nVIDIA GPU with CUDA installed, or a supported AMD GPU with ROCm and ROCm Tensorflow installed.
- clue (available on CRAN)
- parallel (included with R by default)


Package suggested when use this package:

- SIMLR 
Due to some users may [meet fatal errors when using the official SIMLR package](https://github.com/BatzoglouLabSU/SIMLR/issues/47), we built an alternative that remove the optional T-SNE step to  get rid of this issue. To install the unofficial alternative, try `install_github("yulijia/SIMLR", ref = 'master')`.
- SingleCellExperiment
The demo dataset is a SingleCellExperiment object.


The functions in this package are described below.


## encode

**Description**

Generates an encoded subspace of a single-cell RNA-seq expression matrix.

**Usage**

```
encode(dat, seed = 1, max_random_projection = 2048, encoded_dim = 16, hidden_dims = c(128), 
  learning_rate = 0.001, batch_size = 32, epochs = 100, verbose = 1, scale = FALSE,
  genes_as_rows = FALSE)
```

**Arguments**

```
dat                     A matrix, data frame or tibble containing scRNA-seq expression values. By default,
                        genes are assumed to be represented by columns and samples are assumed to be 
                        represented by rows (but see the argument genes_as_rows). NA values are not  
                        supported, but may be replaced by 0s.

seed                    Random seed for initial gene sampling. Currently a seed cannot 
                        be set to reproducibly determine the behaviour of the autoencoder artificial 
                        neural network. 

max_random_projection   Determines the maximum number of genes to be initially sampled prior to 
                        autoencoder training. In practice the number of genes sampled is equal to this 
                        number or 80% of the genes present in the matrix (rounded up), whichever is 
                        smaller.

encoded_dim             The number of features in the encoded data.

hidden_dims             A vector of 1 or more integers, representing the number of nodes in each 
                        successive hidden layer of the encoder half of the autoencoder. Hidden layers in 
                        the decoder use these widths in reverse.

learning_rate           Learning rate for training the artificial neural network.

batch_size              Number of samples per training batch.

epochs                  Number of training epochs.

verbose                 Determines the verbosity of the keras training function. 
                        0: Silent.
                        1: Progress bar.
                        2: One line per epoch.

scale                   If TRUE, gene values are rescaled to a mean of 0 and a standard deviation of 1.

genes_as_rows           If TRUE, rows in the expression matrix are assumed to represent genes and columns 
                        are assumed to represent cells.
```

**Details**

This function accepts a single scRNA-seq expression matrix, randomly samples a number of genes without replacement and trains an autoencoder artificial neural network on the resulting data. The function uses part of this network to encode cell data within a lower-dimensional latent space and returns the encoded matrix. This function does not need to be called directly by the user for clustering (see ensemble_cluster function below), but is provided for greater flexibility.

It is not recommended to run this function in parallel as model training makes use of resources in parallel (CPU cores or GPU, depending on computer setup).

**Value**

An encoded expression matrix wherein cells are represented by rows and latent features are represented by columns.



## ensemble_cluster

**Description**

Generates an ensemble clustering of a single-cell RNA-seq expression matrix.

**Usage**

```
ensemble_cluster(dat, seed = 1, cluster_func = function(x) kmeans(x, centers=5), 
  ensemble_sizes = c(1, 5, 10, 20, 50), cores = 1, ...)
```

**Arguments**

```
dat                 A matrix, data frame or tibble containing scRNA-seq expression values. By default, 
                    genes are assumed to be represented by columns and samples are assumed to be 
                    represented by rows (but see the argument genes_as_rows under the encode function). 
                    NA values are not supported, but may be replaced by 0s.

seed                Used to generate random seeds for the encode function and acts as a random seed 
                    for stochastic clustering functions.

cluster_func        Any function which will accept a matrix (rows as samples, columns as features) and 
                    return a clustering object such as that returned by the kmeans function.

enzemble_sizes      A vector of integers. Number of individual clusterings to be used in each ensemble 
                    clustering returned.

cores               Number of CPU cores to be used in parallel for individual and ensemble clustering.

...                 Optional arguments to be passed to the encode function.
```

**Details**

This function accepts a single scRNA-seq expression matrix. The encode function is used to produce multiple encodings of the data. These are separately clustered using a clustering function optionally provided by the user and produces a set of consensus clusters from these individual clusterings using the clue package, which are returned to the user.

**Value**

A list of length len(ensemble_sizes) containing vectors of consensus clusters per cell. Each ensemble clustering is generated using a number of individual clusterings given by the corresponding element in the ensemble_sizes argument.


## prefilter

A QC function to remove lowly expressed genes 

**Usage**

```
prefilter(table, minReads = 1000, minGene = 100, minCountsperGene = 1, removeZeroGene = T)
```


## estimate_k

Estimated number of clusters via clustering stability metric


**Usage**

```
estimate_k(dat, seed = 1, criteria_method = "NMI", cluster_func = function(x) kmeans(x, centers), krange = 2:25, ensemble_sizes = c(1, 5, 10, 20, 50), cores = 1, ...)
```



**Arguments**

```
dat                 A matrix, data frame or tibble containing scRNA-seq expression values. By default, 
                    genes are assumed to be represented by columns and samples are assumed to be 
                    represented by rows (but see the argument genes_as_rows under the encode function). 
                    NA values are not supported, but may be replaced by 0s.

seed                Used to generate random seeds for the encode function and acts as a random seed 
                    for stochastic clustering functions.

cluster_func        Any function which will accept a matrix (rows as samples, columns as features) and 
                    return a clustering object such as that returned by the kmeans function.

enzemble_sizes      A vector of integers. Number of individual clusterings to be used in each ensemble 
                    clustering returned.

cores               Number of CPU cores to be used in parallel for individual and ensemble clustering.

krange              a range of K that could be evaluated to find the optimal number of clusters

criteria_method     criteria used to measure stability, by default is Normalized Mutual Information (NMI).

...                 Optional arguments to be passed to the encode function.
```



## Example

### scCCESS with K-means

```{r}
library(scCCESS)
library(SingleCellExperiment)

data(sce, package = 'scCCESS')
dat=SingleCellExperiment::counts(sce)

dat=prefilter(dat)
res = estimate_k(dat,
               seed = 1, 
               cluster_func = function(x,centers) { 
                 set.seed(42);
                 kmeans(x, centers)
               },
               criteria_method = "NMI",
               krange = 5:15, ensemble_sizes = 10,
               cores = 8
)

cluster = ensemble_cluster(dat, 
                          seed = 1, 
                          cluster_func = function(x) {
                            set.seed(1)
                            kmeans(x, centers = res$ngroups)
                          }, 
                          cores = 8, 
                          genes_as_rows = T, 
                          ensemble_sizes = 10, 
                          verbose = 0, 
                          scale = F, 
                          batch_size = 64
)


```



### scCCESS with SIMLR

we recommend user to use `SIMLR_Large_Scale`, in which is faster than `SIMLR`.

```{r}
library(scCCESS)
library(SIMLR)

data(sce, package = 'scCCESS')
dat=SingleCellExperiment::counts(sce)

res = estimate_k(dat,
               seed = 1, 
               cluster_func = function(x,centers) {
                 set.seed(42);
                 SIMLR_Large_Scale(t(x), c=centers,kk=15)
               },
               criteria_method = "NMI",
               krange = 5:15, ensemble_sizes = 10,
               cores = 8
)


cluster = ensemble_cluster(dat, 
                          seed = 1, 
                          cluster_func = function(x) {
                            set.seed(1)
                            SIMLR_Large_Scale(t(x), c=res$ngroups,kk=15)
                          }, 
                          cores = 8, 
                          genes_as_rows = T, 
                          ensemble_sizes = 10, 
                          verbose = 0, 
                          scale = F, 
                          batch_size = 64
)
```

## References
Yu, L., Cao, Y., Yang, J. & Yang, P.† (2022) Benchmarking clustering algorithms on estimating the number of cell types from single-cell RNA-sequencing data. _Genome Biology_, 23, 49. [[Full Text](https://doi.org/10.1186/s13059-022-02622-0)]


Geddes, T., Kim, T., Nan, L., Burchfield, J., Yang, J., Tao, D. & Yang, P.† (2019) Autoencoder-based cluster ensembles for single-cell RNA-seq data analysis. BMC Bioinformatics, 20, 660. [[Full Text](https://doi.org/10.1186/s12859-019-3179-5)]
