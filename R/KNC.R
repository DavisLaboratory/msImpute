#' k-nearest class means (KNC)
#'
#' The fraction of k-nearest class means in the original data that are preserved as k-nearest class means in imputed data. KNC 
#' quantifies preservation of the mesoscopic structure after imputation.
#'
#' @param xorigin numeric matrix. The original data. Can contain missing values.
#' @param ximputed numeric matrix. The imputed data.
#' @param class factor. A vector of length number of columns (samples) in the data specifying the class of each sample.
#' @param k  number of nearest class means.
#'
#' @return numeric  The proportion of preserved k-nearest class means in imputed data.
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