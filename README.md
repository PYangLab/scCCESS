# scCCESS

single-cell Consensus Clusters of Encoded Subspaces

This repository stores the R implementation of scCCESS with stability-based approach for estimating number of cell types.


## installation

```
devtools::install_github('scCCESS')
```

Package required to use this package:

- keras (available on CRAN) See https://keras.rstudio.com/reference/install_keras.html for more details - this package can be sped up substantially if the user has a supported nVIDIA GPU with CUDA installed, or a supported AMD GPU with ROCm and ROCm Tensorflow installed.
- clue (available on CRAN)
- parallel (included with R by default)


Package suggested when use this package:

- SIMLR
- SingleCellExperiment


The functions in this package are described below.


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

A pre-filter function to remove lowly expressed genes 

**Usage**

```
prefilter(table, minReads = 1000, minGene = 100, minCountsperGene = 1, removeZeroGene = T)
```


## estimate_k

Estimated number of clusters via clustering stability metric


**Usage**

```
estimate_k(dat, seed = 1, criteria_method = "NMI", cluster_func = function(x) kmeans(x, centers=5), krange = 2:25, ensemble_sizes = c(1, 5, 10, 20, 50), cores = 1, ...)
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

sce=readRDS("demo.rds")

dat=SingleCellExperiment::counts(sce)

dat=prefilter(dat)
ngroups = estimate_k(dat,
                     seed = 1, 
                     cluster_func = function(x,centers) { set.seed(42);kmeans(x, centers)},
                     criteria_method = "NMI",
                     krange = 5:15, ensemble_sizes = 10,
                     cores = 8
)
cluster = ensemble_cluster(dat, seed = 1, cluster_func = function(x) {
  set.seed(1)
  kmeans(x, centers = ngroups)
}, cores = 8, genes_as_rows = T, ensemble_sizes = 10, verbose = 0, 
scale = F, batch_size = 64)

```



### scCCESS with SIMLR


```{r}


library(SIMLR)


ngroups = estimate_k(dat,
                     seed = 1, 
                     cluster_func = function(x,centers) {
                       set.seed(42);
                       SIMLR_Large_Scale(t(x), c=centers,kk=15)
                     },
                     criteria_method = "NMI",
                     krange = 5:15, ensemble_sizes = 10,
                     cores = 8
)


cluster = ensemble_cluster(dat, seed = 1, cluster_func = function(x) {
  set.seed(1)
  SIMLR_Large_Scale(t(x), c=ngroups,kk=15)
}, cores = 8, genes_as_rows = T, ensemble_sizes = 10, verbose = 0, 
scale = F, batch_size = 64)


```
