#' CPD
#'
#' Spearman correlation between pairwise distances in the original data and imputed data. 
#' CPD quantifies preservation of the global structure after imputation.
#'
#' @param xorigin numeric matrix. The original data. Can not contain missing values.
#' @param ximputed numeric matrix. The imputed data. Can not contain missing values.
#'
#' @return numeric 
#'
#' @export
CPD <- function(xorigin, ximputed){
  return(cor(x=as.numeric(dist(t(xorigin))), 
             y = as.numeric(dist(t(ximputed))),
             method = "spearman"))
}