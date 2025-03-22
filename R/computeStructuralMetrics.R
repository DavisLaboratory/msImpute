#' Metrics for the assessment of post-imputation structural preservation 
#'
#' DEPRECATED. For an imputed dataset, it computes within phenotype/experimental condition similarity
#' (i.e. preservation of local structures), between phenotype distances
#' (preservation of global structures), and the Gromov-Wasserstein (GW)
#' distance between original (source) and imputed data.
#'
#' @param x numeric matrix. An imputed data matrix of log-intensity.
#' @param group factor. A vector of biological groups, experimental conditions or
#' phenotypes (e.g. control, treatment).
#' @param y numeric matrix. The source data (i.e. the original log-intensity matrix),
#' preferably subsetted on highly variable peptides (see \code{findVariableFeatures}).
#' @param k numeric. Number of Principal Components used to compute the GW distance.
#' default to 2.
#'
#' @details For each group of experimental conditions (e.g. treatment and control), the group centroid is
#' calculated as the average of observed peptide intensities. Withinness for each group is computed as
#' sum of the squared distances between samples in that group and
#' the group centroid. Betweenness is computed as sum of the squared distances between group centroids.
#' When comparing imputation approaches, the optimal imputation strategy should minimize the within
#' group distances, hence smaller withinness, and maximizes between group distances, hence larger betweenness.
#' The GW metric considers preservation of both local and global structures simultaneously. A small GW distance
#' suggests that imputation has introduced small distortions to global and local structures overall, whereas a
#' large distance implies significant distortions. When comparing two or more imputation methods, the optimal
#' method is the method with smallest GW distance. The GW distance is computed on Principal Components (PCs)
#' of the source and imputed data, instead of peptides. Principal components capture the geometry of the data,
#' hence GW computed on PCs is a better measure of preservation of local and global structures. The PCs in the
#' source data are recommended to be computed on peptides with high biological variance. Hence, users are
#' recommended to subset the source data only on highly variable peptides (hvp) (see \code{findVariableFeatures}).
#' Since the hvp peptides have high biological variance, they are likely to have enough information to discriminate
#' samples from different experimental groups. Hence, PCs computed on those peptides should be representative
#' of the original source data with missing values. If the samples cluster by experimental group in the first
#' couple of PCs, then a choice of k=2 is reasonable. If the desired separation/clustering of samples
#' occurs in later PCs (i.e. the first few PCs are dominated by batches or unwanted variability), then
#' it is recommended to use a larger number of PCs to compute the GW metric.
#' If you are interested in how well the imputed data represent the original data in all possible dimensions,
#' then set k to the number of samples in the data (i.e. the number of columns in the intensity matrix).
#' GW distance estimation requires \code{python}. See example. All metrics are on log scale.
#'
#'
#' @return list of three metrics: withinness (sum of squared distances within a phenotype group),
#' betweenness (sum of squared distances between the phenotypes), and gromov-wasserstein distance (if \code{xna} is not NULL).
#' if \code{group} is NULL only the GW distance is returned. All metrics are on log scale.
#'
#' @references
#' Hediyeh-zadeh, S., Webb, A. I., & Davis, M. J. (2020). MSImpute: Imputation of label-free mass spectrometry peptides by low-rank approximation. bioRxiv.
#'
#' @examples
#' data(pxd010943)
#' y <- log2(data.matrix(pxd010943))
#' y <- y[complete.cases(y),]
#' group <- as.factor(gsub("_[1234]", "", colnames(y)))
#' computeStructuralMetrics(y, group, y=NULL)
#'
#' 
computeStructuralMetrics <- function(x, group=NULL, y = NULL, k=2){
  if(!is.null(group)){
    out <- list(withinness = log(withinness(x, group)),
                betweenness = log(betweenness(x,group)))
 }
  if(!is.null(y)){
    GW <- gromov_wasserstein(x, y, k=k)
    out[['gw_dist']] <- GW[[2]]$gw_dist
  }
  return(out)
}



#' @keywords internal
withinness <- function(x, class_label){
  within_class_dist <- list()
  for(class in class_label){
    centroid <- colMeans(t(x[,class_label==class]))
    within_class_dist[class] <- sum(as.matrix(pdist::pdist(t(x[,class_label==class]), centroid))^2)
  }
  return(unlist(within_class_dist))
}


#' @importFrom stats dist aggregate
#' @keywords internal
betweenness <- function(x, class_label){
  centroids <- aggregate(t(x), list(as.factor(class_label)), mean)
  # the fist column is the group and should be dropped for distance calculation
  return(sum(dist(centroids[,-1])^2))
  #return(lsa::cosine(centroids))

}


#' @importFrom stats prcomp
#' @keywords internal
gromov_wasserstein <- function(x, y, k, min.mean = 0.1){
  if (k > ncol(x)) stop("Number of Principal Components cannot be greater than number of columns (samples) in the data.")
  if (any(!is.finite(x))) stop("Non-finite values (NA, Inf, NaN) encountered in imputed data")
  if (any(!is.finite(y))) stop("Non-finite values (NA, Inf, NaN) encountered in source data")

  means <- rowMeans(x)
  vars <- matrixStats::rowSds(x)

  # Filtering out zero-variance and low-abundance peptides
  is.okay <- !is.na(vars) & vars > 1e-8 & means >= min.mean

  xt <- t(x)
  yt <- t(y)

  # compute PCA
  xt_pca <- prcomp(xt[,is.okay], scale. = TRUE, center = TRUE)
  yt_pca <- prcomp(yt, scale. = TRUE, center = TRUE)

  C1 <- yt_pca$x[, seq_len(k)]
  C2 <- xt_pca$x[, seq_len(k)]


  cat("Computing GW distance using k=", k, "Principal Components\n")
  # reticulate::source_python(system.file("python", "gw.py", package = "msImpute"))
  # return(gw(C1,C2, ncol(x)))
}




