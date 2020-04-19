#' Select features with high dropout rate
#'
#' Fits a linear model to peptide dropout rate against peptide abundance. The selected features (peptides) can be
#' used to determine if data is Missing Not At Random (MNAR). Users should note that \code{msImpute} assumes peptides
#' are Missing At Random (MAR).
#'
#' @param object Numeric matrix or \code{MAList} object from \link{limma} where missing values are denoted by NA.
#' Rows are peptides, columns are samples.
#' @param n_features Numeric, number of features with high dropout rate. 500 by default.
#' @param suppress_plot Logical show plot of dropouts vs abundances.
#'
#' @return A data frame with a logical column denoting the selected features
#'
#' @examples
#' set.seed(101)
#' n=200
#' p=100
#' J=50
#' np=n*p
#' missfrac=0.3
#' x=matrix(rnorm(n*J),n,J)%*%matrix(rnorm(J*p),J,p)+matrix(rnorm(np),n,p)/5
#' ix=seq(np)
#' imiss=sample(ix,np*missfrac,replace=FALSE)
#' xna=x
#' xna[imiss]=NA
#' selectFeatures(xna, n_features=500,  suppress_plot=FALSE)
#' @export
selectFeatures <- function(object, n_features=500, suppress_plot = FALSE) {

  if(is(object,"MAList")){
    x <- object$E
  }else{
    x <- object
  }

  if(is.null(rownames(x))) stop("No row names in input. Please provide input with named rows.")
  AveExpr <- rowMeans(x, na.rm = TRUE)
  dropout <- rowMeans(is.na(x))

  linear_fit <- lm(dropout ~ AveExpr)
  resids <- residuals(linear_fit)
  lin_res_o <- order(resids, decreasing = TRUE)

  cols <- rep("#3E71A8", length(resids))
  cols[lin_res_o[1:n_features]] <- "#DE1A1A"

  if(!suppress_plot){
    plot(x = AveExpr, y = dropout, pch = 16,
         cex = 0.5, col = cols, main = paste("Top ",n_features," high droupout peptides", sep =""))
    abline(linear_fit)
  }

  hdrp <- data.frame(name = rownames(x), AveExpr = AveExpr, dropout = dropout ,
                         residual = resids, msImpute_feature=FALSE)
  hdrp$msImpute_feature[lin_res_o[1:n_features]] <- TRUE
  hdrp <- data.table::as.data.table(hdrp)

  return(hdrp)

}
