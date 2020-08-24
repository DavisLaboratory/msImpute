#' k-nearest neighbour (KNN)
#'
#' The fraction of k-nearest neighbours in the original data that are preserved as k-nearest neighbours in imputed data.
#' KNN quantifies preservation of the local, or microscopic structure.
#' Requires complete datasets - for developers/use in benchmark studies only.
#'
#' @param xorigin numeric matrix. The original log-intensity data. Can not contain missing values.
#' @param ximputed numeric matrix. The imputed log-intensity data. Can not contain missing values.
#' @param k  number of nearest neighbours. default to k=3.
#'
#' @return numeric  The proportion of preserved k-nearest neighbours in imputed data.
#' @examples
#' data(pxd007959)
#' y <- pxd007959$y
#' y <- y[complete.cases(y),]
#' # for demonstration we use same y for xorigin and ximputed
#' KNN(y, y)
#'
#'
#' @export
KNN <- function(xorigin, ximputed, k=3){

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
