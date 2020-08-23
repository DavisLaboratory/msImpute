#' Standardize a matrix to have optionally row means zero and variances one, and/or column means zero and variances one.
#'
#'
#' @param object numeric matrix giving log-intensity where missing values are denoted by NA. Rows are peptides, columns are samples.
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
#' Data is assumed to be normalised and log-transformed.
#'
#' @return
#' A list of two components: E and E.scaled. E contains the input matrix, E.scaled contains the scaled data
#'
#'
#' @examples
#' set.seed(101)
#' n=12000
#' p=10
#' J=5
#' np=n*p
#' missfrac=0.3
#' x=matrix(rnorm(n*J,mean = 5,sd = 0.2),n,J)%*%matrix(rnorm(J*p, mean = 5,sd = 0.2),J,p)+
#'   matrix(rnorm(np,mean = 5,sd = 0.2),n,p)/5
#' ix=seq(np)
#' imiss=sample(ix,np*missfrac,replace=FALSE)
#' xna=x
#' xna[imiss]=NA
#' keep <- (rowSums(!is.na(xna)) >= 4)
#' xna <- xna[keep,]
#' xna <- scaleData(xna)
#' @seealso selectFeatures, msImpute
#' @references
#' Hastie, T., Mazumder, R., Lee, J. D., & Zadeh, R. (2015). Matrix completion and low-rank SVD via fast alternating least squares. The Journal of Machine Learning Research, 16(1), 3367-3402.
#' @references
#' Hediyeh-zadeh, S., Webb, A. I., & Davis, M. J. (2020). MSImpute: Imputation of label-free mass spectrometry peptides by low-rank approximation. bioRxiv.
#' @importFrom methods is
#' @export
scaleData <- function(object, maxit = 20, thresh = 1e-09, row.center = TRUE, row.scale =TRUE,
                      col.center = TRUE, col.scale = TRUE, trace = FALSE){
  if(is(object,"MAList")){
    x <- object$E
  }else{
    x <- object
  }
  if(any(is.nan(x) | is.infinite(x))) stop("Inf or NaN values encountered.")
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
