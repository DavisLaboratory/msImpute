#' Plot mean-CV^2 trend
#'
#' For each peptide, the squares of coefficient of variations are computed and plotted against average log-intensity.
#' Additionally, a loess trend is fitted to the plotted values.
#' Outlier observations (possibly originated from incorrect match between runs), are detected and highlighted.
#' Users can use this plot as a diagnostic plot to determine if filtering by average intensity is required.
#'
#' @details
#' Outliers are determined by computing the RBF kernels, which reflect the chance that an observed point
#' belong to the dataset (i.e. is close enough in distance to other data points). Users can determine the cut-off
#' for intensity-based filtering with respect to the mean log-intensity of the outlier points.
#'
#' @param y numeric matrix of log-intensity
#' @param trend logical. Should a loess trend be fitted to CV^2 and mean values. Default to TRUE.
#' @param main character string. Title of the plot. Default  to NULL
#'
#' @return A plot is created on the current graphics device.
#' @examples
#' data(pxd010943)
#' y <- pxd010943
#' y <- log2(y)
#' ppCV2 <- plotCV2(y)
#'
#' @importFrom limma loessFit
#' @importFrom matrixStats rowSds
#' @importFrom graphics plot lines points
#' @export
plotCV2 <- function(y, trend = TRUE, main=NULL){
  A <- rowMeans(y, na.rm = TRUE)
  CV <- (matrixStats::rowSds(data.matrix(y), na.rm = TRUE)/A)^2
  res <- data.frame(mean = A, CV = CV)
  plot(A, CV, cex = 0.3, pch = 16,
       xlab="Average log-intensity", ylab=expression("CV"^2), main=main)
  if(trend){
    fit <- limma::loessFit(CV, A)
    o <- order(A)
    lines(A[o], fit$fitted[o], lwd =2, col = "red")
  }

  return(res)
}
