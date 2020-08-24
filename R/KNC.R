#' k-nearest class means (KNC)
#'
#' The fraction of k-nearest class means in the original data that are preserved as k-nearest class means in imputed data. KNC
#' quantifies preservation of the mesoscopic structure after imputation.
#' Requires complete datasets - for developers/use in benchmark studies only.
#'
#' @param xorigin numeric matrix. The original log-intensity data. Can contain missing values.
#' @param ximputed numeric matrix. The imputed log-intensity data.
#' @param class factor. A vector of length number of columns (samples) in the data specifying the class/label (i.e. experimental group) of each sample.
#' @param k  number of nearest class means. default to k=3.
#'
#' @return numeric  The proportion of preserved k-nearest class means in imputed data.
#'
#' @examples
#' data(pxd007959)
#' y <- pxd007959$y
#' y <- y[complete.cases(y),]
#' # for demonstration we use same y for xorigin and ximputed
#' KNC(y, y, class = as.factor(pxd007959$samples$group))
#'
#' @export
KNC <- function(xorigin, ximputed, class, k=3){
  class_means_org <- list()
  for(G in unique(class)){
    class_means_org[[G]] <- rowMeans(xorigin[,class ==G], na.rm = TRUE)
  }
  NN_org <- FNN::get.knn(t(data.frame(class_means_org)), k = k)
  KNC_org <- NN_org$nn.index

  class_means_amp <- list()
  for(G in unique(class)){
    class_means_amp[[G]] <- rowMeans(ximputed[,class==G])
  }

  NN_amp <- FNN::get.knn(t(data.frame(class_means_amp)), k = k)
  KNC_amp <- NN_amp$nn.index
  pmeans <- c()
  for(i in 1:length(levels(class))){
    pmeans <- c(pmeans, mean(KNC_amp[i,] %in% KNC_org[i,]))
  }
  return(mean(pmeans))
}
