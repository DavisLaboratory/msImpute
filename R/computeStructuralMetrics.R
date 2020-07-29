#' Metrics for the assessment of post-imputation structural preservation
#'
#' For an imputed dataset, it computes within phenotype/experimental condition similarity (i.e. preservation of local structures)
#' and between phenotype distances (preservation of global structures). Metrics are on log scale.
#'
#' @param x numeric matrix. An imputed data matrix.
#' @param group factor. A vector of biological groups, experimental conditions or phenotypes (e.g. control, treatment).
#'
#' @details For each group of experimental conditions (e.g. treatment and control), the group centroid is calculated as the average
#' of observations. Withinness for each group is computed as sum of the squared distances between samples in that group and the group centroid. Betweenness is computed as
#' sum of the squared distances between group centroids. When comparing imputation approaches, the optimal imputation strategy should minimize the within
#' group distances (hence smaller withinness), and maximizes between group distances, hence larger betweenness.
#'
#' @return list of two metrics: withinness (sum of squared distances within a phenotype group on log scale),
#' betweenness (sum of squared distances between the phenotypes on log scale)
#'
#' @export
computeStructuralMetrics <- function(x, group){
  return(list(withinness = log(withinness(x, group)),
              betweenness = log(betweenness(x,group))))

}



#' @export
withinness <- function(x, class_label){
  within_class_dist <- list()
  for(class in class_label){
  centroid <- colMeans(t(x[,class_label==class]))
  within_class_dist[class] <- sum(as.matrix(pdist::pdist(t(x[,class_label==class]), centroid))^2)
  }
  return(unlist(within_class_dist))
}

#' @export
betweenness <- function(x, class_label){
  centroids <- sapply(unique(class_label), FUN=function(class){
    colMeans(t(x[,class_label==class]))
  })

  return(sum(dist(centroids)^2))
  #return(lsa::cosine(centroids))

}
