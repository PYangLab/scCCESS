#' prefilter
#'
#' Quality control: to filter the data to only include true cells that are of high quality
#'
#'
#' @param table input scRNAseq matrix, where each row represents a gene
#' and each column represents a single cell with a raw count for every row (gene) in the matrix.
#' @param minReads minimum reads count (UMI) per cell
#' @param minGene minimum gene that have counts (UMI) more than \code{minCountsperGene}
#' @param minCountsperGene minimum counts (UMI) are required for each gene
#' @param removeZeroGene remove genes that has zero counts were detected in the cell
#'
#' @return a filtered scRNAseq expression matrix
#'
#'@examples
#'
#' # library(SingleCellExperiment)
#' # library(scCCESS)
#' # data("sce", package = "scCCESS")
#' # dat=SingleCellExperiment::counts(sce)
#' # dat.filtered=prefilter(dat)
#'
#' @importFrom Matrix colSums rowSums
#'
#' @export
#'
prefilter = function(table, minReads = 1000, minGene = 100, minCountsperGene = 1, removeZeroGene = T){
  nGene = colSums(table >= minCountsperGene)
  nReads = colSums(table)
  keepLibrary = nReads >= minReads
  keepGene = nGene >= minGene
  ZeroGene = rowSums(table)
  if (removeZeroGene == TRUE ){
    removeGene = ZeroGene = 0
    dat.filtered = table[!removeGene, keepLibrary & keepGene]
    dat.filtered = log2( dat.filtered + 1 )
  }else{
    dat.filtered = table[, keepLibrary & keepGene]
    dat.filtered = log2( dat.filtered + 1 )
  }
  cat("total removed genes ", sum(removeGene), "\t; total keeped cells ", sum(keepLibrary & keepGene), "\n")
  return(dat.filtered)
}

#' Generates an encoded subspace of a single-cell RNA-seq expression matrix
#'
#' This function accepts a single scRNA-seq expression matrix,
#' randomly samples a number of genes without replacement
#' and trains an autoencoder artificial neural network on the resulting data.
#' The function uses part of this network to encode cell data within
#' a lower-dimensional latent space and returns the encoded matrix.
#' This function does not need to be called directly by the user for clustering (see ensemble_cluster function below),
#' but is provided for greater flexibility.
#' It is not recommended to run this function in parallel
#' as model training makes use of resources in parallel (CPU cores or GPU, depending on computer setup).
#'
#'
#' @param dat  A raw data matrix, data frame or tibble containing scRNA-seq expression values. By default,
#' genes are assumed to be represented by columns and samples are assumed to be
#' represented by rows (but see the argument genes_as_rows). NA values are not
#' supported, but may be replaced by 0s. It could be a filtered matrix by prefilter() funtion.
#' @param seed Random seed for initial gene sampling. Currently a seed cannot
#' be set to reproducibly determine the behaviour of the autoencoder artificial
#' neural network.
#' @param max_random_projection Determines the maximum number of genes to be initially sampled prior to
#' autoencoder training. In practice the number of genes sampled is equal to this
#' number or 80% of the genes present in the matrix (rounded up), whichever is
#' smaller.
#' @param encoded_dim  The number of features in the encoded data.
#' @param hidden_dims A vector of 1 or more integers, representing the number of nodes in each
#' successive hidden layer of the encoder half of the autoencoder. Hidden layers in
#' the decoder use these widths in reverse.
#' @param learning_rate Learning rate for training the artificial neural network.
#' @param batch_size Number of samples per training batch.
#' @param epochs Number of training epochs.
#' @param verbose Determines the verbosity of the keras training function, 0 is silent, 1 is an animated progress bar, 2 is message, default is 2
#' @param scale If TRUE, gene values are rescaled to a mean of 0 and a standard deviation of 1.
#' @param genes_as_rows If TRUE, rows in the expression matrix are assumed to represent genes and columns
#' assumed to represent cells.
#' @param cores number of cores
#'
#' @return An encoded expression matrix wherein cells are represented by rows and latent features are represented by columns.
#'
#'
#' @examples
#'
#' # library(SingleCellExperiment)
#' # library(scCCESS)
#' # data("sce", package = "scCCESS")
#' # dat=SingleCellExperiment::counts(sce)
#' # dat.filtered=prefilter(dat)
#' # encoding = encode(dat.filtered, seed = 1,
#' #                   genes_as_rows = T, ensemble_sizes = 10,
#' #                   verbose = 0, scale = F, batch_size = 64,
#' #                   cores=8)
#'
#' @import keras tensorflow
#' @importFrom methods is
#' @importFrom stats kmeans median predict sd
#' @importFrom utils packageVersion
#'
#' @export
encode = function(dat, seed = 1, max_random_projection = 2048, encoded_dim = 16, hidden_dims = c(128), learning_rate = 0.001, batch_size = 32, epochs = 100, verbose = 2, scale = FALSE, genes_as_rows = FALSE,cores=1) {
  if (verbose[1] %in% 0:2) {
    verbose = verbose[1]
  } else {
    verbose = 1
  }

  set.seed(seed)
  was_data_frame=is.data.frame(dat)
  if (is.data.frame(dat)||is(dat,"Matrix")) dat = as.matrix(dat)
  if (!is(dat,"matrix")) stop("Input data must be dataframe or matrix.")

  # Transpose
  if (genes_as_rows) dat = t(dat)

  # Strip row and column names
  datrows = rownames(dat)
  rownames(dat) = NULL
  colnames(dat) = NULL

  # Scale columns
  if (scale) {
    dat = apply(dat, 2, function(x) (x - mean(x)) / sd(x))
  }

  # Perform random projection
  num_input_features = ncol(dat)
  final_proj_dim = min(max_random_projection, ceiling(0.8 * num_input_features))
  random_proj_cols = sample(num_input_features, size = final_proj_dim, replace = FALSE)
  dat = dat[, random_proj_cols]

  # Clear deep learning graph
  k_clear_session()

  # Construct encoder network
  tns = encoder_input = layer_input(shape = final_proj_dim)
  for (layer_width in hidden_dims) {
    tns = layer_dense(tns, units = layer_width)
    tns = layer_activation_leaky_relu(tns, alpha = 0.01)
  }
  tns = layer_dense(tns, units = encoded_dim)
  encoder = keras_model(inputs = encoder_input, outputs = tns)

  # Construct decoder network
  tns = decoder_input = layer_input(shape = encoded_dim)

  for (layer_width in rev(hidden_dims)) {
    tns = layer_dense(tns, units = layer_width)
    tns = layer_activation_leaky_relu(tns, alpha = 0.01)
  }

  tns = layer_dense(tns, units = final_proj_dim)
  decoder = keras_model(inputs = decoder_input, outputs = tns)

  # Combine networks
  tns = ae_input = layer_input(final_proj_dim)
  tns = decoder(encoder(tns))
  autoencoder = keras_model(inputs = ae_input, outputs = tns)
  if(packageVersion("keras")>="2.6.0"){
    compile(autoencoder, optimizer = optimizer_adam(learning_rate = learning_rate), loss = 'mean_squared_error')
  }else{
    compile(autoencoder, optimizer = optimizer_adam(lr = learning_rate), loss = 'mean_squared_error')
  }


  # Fit autoencoder model
  fit(autoencoder, dat, dat, batch_size = batch_size, epochs = epochs, verbose = verbose)

  # Encode input data, return rownames and return as original data type
  reduced_data = predict(encoder, dat, batch_size = batch_size)
  rownames(reduced_data) = datrows
  # if (genes_as_rows) reduced_data %<>% t
  if (was_data_frame) reduced_data = as.data.frame(reduced_data)

  return(reduced_data)
}


#' Generates an ensemble clustering of a single-cell RNA-seq expression matrix
#'
#'
#'
#' @param dat  A matrix, data frame or tibble containing scRNA-seq expression values. By default,
#' genes are assumed to be represented by columns and samples are assumed to be
#' represented by rows (but see the argument genes_as_rows under the encode function).
#' NA values are not supported, but may be replaced by 0s.
#' @param seed  Used to generate random seeds for the encode function and acts as a random seed
#' for stochastic clustering functions.
#' @param cluster_func Any clustering function which will accept a matrix (rows as samples, columns as features).
#' It could be any function that need provide a range of k, such as K-means, SIMLR, C-means, K-medoids, etc.
#' @param ensemble_sizes A vector of integers. Number of individual clusterings to be used in
#' each ensemble clustering returned.
#' @param cores  Number of CPU cores to be used in parallel for individual and ensemble clustering.
#' @param ...   Optional arguments to be passed to the encode function.
#'
#' @return A list of length len(ensemble_sizes) containing vectors of consensus clusters per cell.
#' Each ensemble clustering is generated using a number of individual clusterings
#' given by the corresponding element in the ensemble_sizes argument.
#'
#' @examples
#'
#' # library(SingleCellExperiment)
#' # library(scCCESS)
#' # data("sce", package = "scCCESS")
#' # dat=SingleCellExperiment::counts(sce)
#' # cluster = ensemble_cluster(dat, seed = 1, cluster_func = function(x) { set.seed(42)
#' #                            kmeans(x, centers = 7)},
#' #                            cores = 8, genes_as_rows = T,
#' #                            ensemble_sizes = 10, verbose = 0,
#' #                            scale = F, batch_size = 64)
#'
#' @export
#'
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @importFrom stats kmeans median predict sd
#' @import clue
#'
ensemble_cluster = function(dat, seed = 1, cluster_func = function(x) kmeans(x, centers=5), ensemble_sizes = c(1, 5, 10, 20, 50), cores = 1, ...) {

  individual_encodings = lapply(1:max(ensemble_sizes), function(i) {
    ae_seed = (seed - 1) * max(ensemble_sizes) + i
    cat("Seed ", ae_seed, "\tEnsemble ", i, "\n")
    encoding = encode(dat, seed = ae_seed,cores=cores, ...)
    cat("\n")
    return(encoding)
  })

  cat("Generating individual clusters...\n")
  individual_clusters = mclapply(individual_encodings, cluster_func, mc.cores = cores)

  cat("Generating consensus clusters...\n")
  ensemble_clusters = mclapply(ensemble_sizes, function(size) {
    if (size == 1) {
      if(grepl("SIMLR",deparse1(cluster_func))){
       return(sample(individual_clusters, 1)[[1]]$y$cluster)
      }else{
        return(sample(individual_clusters, 1)[[1]]$cluster)
      }
    } else {
      consensus = sample(individual_clusters, size, replace = FALSE)
      if(grepl("SIMLR",deparse1(cluster_func))){
        consensus=lapply(consensus, function(x)return(x$y))
      }
      consensus = cl_consensus(consensus, method = "HE")

      consensus = apply(as.matrix(consensus[[1]]), 1, function(row_vals) which(row_vals == 1))
      return(consensus)
    }
  }, mc.cores = cores)

  if(length(ensemble_sizes)==1){
    return(unlist(ensemble_clusters))
  }else{
    names(ensemble_clusters) = as.character(ensemble_sizes)
    return(ensemble_clusters)
  }

}


#' Clustering concordance by ARI, NMI, FM or Jaccard index
#'
#' @param x  clusters A
#' @param y clusters B
#' @param criteria measures include ARI, NMI
#'
#' @return the similarity value of two clusters.
#'
#' @examples
#'
#' # criteria_compare(x=c(1,2,2,2,1,2),y=c(1,2,2,3,1,2),"ARI")
#'
#' @keywords internal
#'
#' @import aricode
#'
criteria_compare = function(x, y, criteria) {
  res = switch(criteria,
               "ARI" = round(ARI(x, y), digits = 5),
               "NMI" = round(NMI(x, y), digits = 5)#,
               # "FM" =  round(clusterCrit::extCriteria(x, y, "Folkes_Mallows")$folkes_mallows, digits = 5),
               # "Jaccard" =  round(clusterCrit::extCriteria(x, y, "Jaccard")$jaccard, digits = 5)
  )

  return(res)
}



#' Estimated number of clusters via stability metrics
#'
#'
#'
#' @param dat A matrix, data frame or tibble containing scRNA-seq expression values. By default,
#' genes are assumed to be represented by columns and samples are assumed to be
#' represented by rows (but see the argument genes_as_rows under the encode function).
#' NA values are not supported, but may be replaced by 0s.
#' @param seed random seed
#' @param criteria_method stability metrics to be used, default is NMI, this should be one of "NMI", "ARI", "Jaccard".
#' @param cluster_func Any clustering function which will accept a matrix (rows as samples, columns as features).
#' It could be any function that need provide a range of k, such as K-means, SIMLR, C-means, K-medoids, etc.
#' @param krange a set of k (e.g. 2:10) that could be estimated as the number of cell types.
#' @param ensemble_sizes A vector of integers. Number of individual clusterings to be used in
#' each ensemble clustering returned.
#' @param cores Number of CPU cores to be used in parallel for individual and ensemble clustering.
#' @param ... Optional arguments to be passed to the encode function.
#'
#' @return the number of cell types, and the pairwise stability score of each k
#'
#' @examples
#'
#' # library(SingleCellExperiment)
#' # library(scCCESS)
#' # data("sce", package = "scCCESS")
#' # dat=SingleCellExperiment::counts(sce)
#' # k = estimate_k(dat,
#' #                seed = 1,
#' #                cluster_func = function(x,centers) { set.seed(42);kmeans(x, centers)},
#' #                criteria_method = "NMI",
#' #                krange = 2:20, ensemble_sizes = 10,
#' #                cores = 8)
#'
#'
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @importFrom stats kmeans median predict sd
#'
#' @export
#'
#'
estimate_k = function(dat, seed = 1, criteria_method = "NMI", cluster_func = function(x) kmeans(x, centers=5), krange = 2:25, ensemble_sizes = c(1, 5, 10, 20, 50), cores = 1, ...){

  individual_encodings = lapply(1:max(ensemble_sizes), function(i) {
    ae_seed = (seed - 1) * max(ensemble_sizes) + i
    cat("Seed ", ae_seed, "\tEnsemble ", i, "\n")
    encoding = encode(dat, seed = ae_seed, genes_as_rows = T, verbose = 0, scale = F,
                      batch_size = 64,cores=cores ,...)
    cat("\n")
    return(encoding)
  })

  mc = lapply(individual_encodings, function(X) mclapply(krange, cluster_func, mc.cores = cores ,x = X))
  df =list()
  criteria.result = c()
  criteria.list=list()
  for (i in 1:length(krange)) {
    df[[i]] = do.call(rbind, lapply(mc, function(x) {
      if(grepl("SIMLR",deparse1(cluster_func))){
        return(x[[i]]$y$cluster)
      }else{
        return(x[[i]]$cluster)
      }

    }))

    l = 1
    criteria = c()
    for (a in 1:(max(ensemble_sizes) - 1)) {
      for (b in (a + 1):max(ensemble_sizes)) {
        criteria[l] = criteria_compare(df[[i]][a, ], df[[i]][b, ],criteria_method)
        l = l + 1
      }
    }
    criteria.list[[i]] = criteria
    criteria.result[i] = median(criteria)
  }


  stability.df=data.frame(do.call(cbind,criteria.list))
  colnames(stability.df)=krange
  ngroups = krange[which.max(criteria.result)]
  return(list(ngroups=ngroups,stability.df=stability.df))
}
