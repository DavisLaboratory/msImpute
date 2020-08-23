#' Find highly variable peptides
#'
#' For each peptide, the total variance is decomposed into biological and technical variance using package \code{scran}
#' @param y numeric matrix giving log-intensity. Can contain NA values.
#'
#' @return A data frame where rows are peptides and columns contain estimates of biological and technical variances. Peptides are ordered by biological variance.
#'
#' @details A loess trend is fitted to total sample variances and mean intensities. For each peptide, the biological variance is then
#' computed by subtracting the estimated technical variance from the loess fit from the total sample variance.
#'
#' @seealso computeStructuralMetrics
#'
#' @export
#' @importFrom scran trendVar decomposeVar
#' @importFrom graphics lines plot
findVariableFeatures <- function(y){
  fit <- trendVar(y)
  results <- decomposeVar(y, fit)
  plot(results$mean, results$total)
  o <- order(results$mean)
  lines(results$mean[o], results$tech[o], col="red", lwd=2)
  results <- as.data.frame(results)
  top.dec <- results[order(results$bio, decreasing=TRUE), ]
  return(top.dec)

}
