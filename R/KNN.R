#' k-nearest neighbour (KNN)
#'
#' The fraction of k-nearest neighbours in the original data that are preserved as k-nearest neighbours in imputed data.
#' KNN quantifies preservation of the local, or microscopic structure.
#'
#' @param xorigin numeric matrix. The original data. Can not contain missing values.
#' @param ximputed numeric matrix. The imputed data. Can not contain missing values.
#' @param k  number of nearest neighbours. default to k=3.
#'
#' @return numeric  The proportion of preserved k-nearest neighbours in imputed data.
#'
#' @export
KNN <- function(xorigin, ximputed, class, k=3){
  
  NN_org <- FNN::get.knn(t(xorigin), k = k)
  KNC_org <- NN_org$nn.index
  
  
  NN_amp <- FNN::get.knn(t(ximputed), k = k)
  KNC_amp <- NN_amp$nn.index
  pmeans <- c()
  for(i in 1:ncol(xorigin)){
    pmeans <- c(pmeans, mean(KNC_amp[i,] %in% KNC_org[i,]))
  }
  return(mean(pmeans))
}