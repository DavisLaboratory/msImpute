#' Standardize a matrix to have optionally row means zero and variances one, and/or column means zero and variances one.
#'
#'
#' @param object numeric matrix where missing values are denoted by NA. Rows are peptides, columns are samples.
#' @param maxit numeric. maximum iteration for the algorithm to converge (default to 20). When both row and column centering/scaling is requested, iteration may be necessary.
#' @param thresh numeric. Convergence threshold (default to 1e-09).
#' @param row.center logical. if row.center==TRUE (the default), row centering will be performed resulting in a matrix with row means zero. If row.center is a vector, it will be used to center the rows. If row.center=FALSE nothing is done.
#' @param row.scale if row.scale==TRUE, the rows are scaled (after possibly centering, to have variance one. Alternatively, if a positive vector is supplied, it is used for row centering.
#' @param col.center Similar to row.center
#' @param col.scale  Similar to row.scale
#' @param trace logical. With trace=TRUE, convergence progress is reported, when iteration is needed.
#'
#' @details
#' Standardizes rows and/or columns of a matrix with missing values, according to the \code{biScale} algorithm in Hastie et al. 2015.
#'
#' @return
#' A list of two components: E and E.scaled. E contains the input matrix, E.scaled contains the scaled data
#'
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
#' xna <- scaleData(xna)
#' @seealso selectFeatures, msImpute
#' @export
scaleData <- function(object, maxit = 20, thresh = 1e-09, row.center = TRUE, row.scale =TRUE,
                      col.center = TRUE, col.scale = TRUE, trace = FALSE){
  if(is(object,"MAList")){
    x <- object$E
  }else{
    x <- object
  }

  if(any(rowSums(!is.na(x)) <= 3)) stop("Peptides with excessive NAs are detected. Please revisit your fitering step. At least 4 non-missing measurements are required for any peptide.")
  if(any(x < 0, na.rm = TRUE)){
    warning("Negative values encountered in imputed data. Please consider revisting the filtering and/or normalisation steps, if appropriate.")
  }

  cat("bi-scaling ...\n")
  xnas <- softImpute::biScale(x, maxit = maxit, thresh = thresh, row.center = row.center, row.scale =row.scale,
                              col.center = col.center, col.scale = col.scale, trace = trace)
  cat("data scaled \n")

  return(list(E = object, E.scaled = xnas))

  # if(is(object,"MAList")) {
  #   object$scaledData <- xnas
  #   return(object)
  # }else{
  #   return(list(object = x, scaledData = xnas))
  # }


}
