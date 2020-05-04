#' Peptide-level imputation in mass spectrometry label-free proteomics by low-rank approximation
#'
#' Returns a completed peptide intensity matrix where missing values (NAs) are imputated
#' by low-rank approximation of the input matrix. Non-NA entries remain unmodified. \code{msImpute} requires at least 4
#' non-missing measurements per peptide across all samples. It is assumed that peptide intensities  (DDA), or MS1/MS2 normalised peak areas (DIA),
#' are log2-transformed and normalised (e.g. quantile normalisation).
#'
#' \code{msImpute} operates on the softImpute-ALS algorithm.
#' For more details on the underlying algorithm, please see \code{\link[softImpute]{softImpute}} package.
#'
#' @param object Numeric matrix or \code{MAList} object from \link{limma} where missing values are denoted by NA. Rows are peptides, columns are samples.
#' @param rank.max Numeric. This restricts the rank of the solution. is set to min(dim(\code{object})-1) by default.
#' @param lambda Numeric. Nuclear-norm regularization parameter. Controls the low-rank property of the solution
#' to the matrix completion problem. By default, it is determined at the scaling step. If set to zero
#' the algorithm reverts to "hardImputation", where the convergence will be slower.
#' @param thresh Numeric. Convergence threshold. Set to 1e-05, by default.
#' @param maxit Numeric. Maximum number of iterations of the algorithm before the algorithm is converged. 100 by default.
#' @param trace.it Logical. Prints traces of progress of the algorithm.
#' @param warm.start List. A SVD object can be used to initialize the algorithm instead of random initialization.
#' @param final.svd  Logical. Shall final SVD object be saved? The solutions to the matrix completion problems are computed from U, D and V components of final SVD.
#' @return Missing values are imputed by low-rank approximation of the input matrix. If input is a numeric matrix,
#' a numeric matrix of identical dimensions is returned. If \code{x} is a \code{MAList} object, the \code{E} component is
#' replaced with the completed matrix, and the updated \code{MAList} object is returned. Non-NA entries remain unmodified.
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
#' msImpute(object=xna)
#' @export
msImpute <- function(object, rank.max = NULL, lambda = NULL, thresh = 1e-05,
                     maxit = 100, trace.it = FALSE, warm.start = NULL, final.svd = TRUE) {

  if(is(object,"MAList")){
    x <- object$E
  }else{
    x <- object
  }

  if(any(rowSums(!is.na(x)) <= 3)) stop("Peptides with excessive NAs are detected. Please revisit your fitering step. At least 4 non-missing measurements are required for any peptide.")
  cat("bi-scaling ...\n")
  xnas <- softImpute::biScale(x)
  cat("data scaled \n")
  if(is.null(rank.max)) rank.max <- min(dim(x) - 1)
  cat("maximum rank is", rank.max, "\n")
  cat("computing lambda0 ... \n")
  if(is.null(lambda)) lambda <- softImpute::lambda0(x)
  cat("lambda0 is", lambda, "\n")
  cat("fit the low-rank model ... \n")
  fit <- softImpute::softImpute(xnas,rank=rank.max,lambda=lambda, type = "als", thresh = thresh,
                  maxit = maxit, trace.it = trace.it, warm.start = warm.start, final.svd = final.svd)
  cat("model fitted. \nImputting missing entries ... \n")
  ximp <- softImpute::complete(x, fit)
  cat("Imputation completed \n")  # need to define a print method for final rank model fitted

  if(any(object < 0)){
    warning("Negative values encountered in imputed data. Please consider revising filtering and/or normalisation steps.")
  }

  if(is(object,"MAList")) {
    object$E <- ximp
    return(object)
  }else{
    return(ximp)
  }


}
