#' Imputation of peptide log-intensity in mass spectrometry label-free proteomics by low-rank approximation
#'
#' Returns a completed matrix of peptide log-intensity where missing values (NAs) are imputated
#' by low-rank approximation of the input matrix. Non-NA entries remain unmodified. \code{msImpute} requires at least 4
#' non-missing measurements per peptide across all samples. It is assumed that peptide intensities  (DDA), or MS1/MS2 normalised peak areas (DIA),
#' are log2-transformed and normalised (e.g. by quantile normalisation).
#'
#' @details
#'
#' \code{msImpute} operates on the \code{softImpute-als} algorithm in \code{\link[softImpute]{softImpute}} package.
#' The algorithm estimates a low-rank matrix ( a smaller matrix
#' than the input matrix) that approximates the data with a reasonable accuracy. \code{SoftImpute-als} determines the optimal
#' rank of the matrix through the \code{lambda} parameter, which it learns from the data.
#' This algorithm is implemented in \code{method="v1"}.
#' In v2 we have used a information theoretic approach to estimate the optimal rank, instead of relying on \code{softImpute-als}
#' defaults. Similarly, we have implemented a new approach to estimate \code{lambda} from the data. Low-rank approximation
#' is a linear reconstruction of the data, and is only appropriate for imputation of MAR data. In order to make the
#' algorithm applicable to MNAR data, we have implemented \code{method="v2-mnar"} which imputes the missing observations
#' as weighted sum of values imputed by msImpute v2 (\code{method="v2"}) and random draws from a Gaussian distribution.
#' Missing values that tend to be missing completely in one or more experimental groups will be weighted more (shrunken) towards
#' imputation by sampling from a Gaussian parameterised by smallest observed values in the sample (similar to minProb, or
#' Perseus). However, if the missing value distribution is even across the samples for a peptide, the imputed values
#' for that peptide are shrunken towards
#' low-rank imputed values. The judgment of distribution of missing values is based on the EBM metric implemented in
#' \code{selectFeatures}, which is also a information theory measure.
#'
#'
#' @param y Numeric matrix giving log-intensity where missing values are denoted by NA. Rows are peptides, columns are samples.
#' @param method character. Allowed values are \code{"v2"} for \code{msImputev2} imputation (enhanced version) for MAR.
#' \code{method="v2-mnar"} (modified low-rank approx for MNAR), and \code{"v1"} initial release of \code{msImpute}
#' @param group character or factor vector of length \code{ncol(y)}
#' @param alpha numeric. The weight parameter. Default to 0.2. Weights the MAR-imputed distribution in the imputation scheme.
#' @param rank.max Numeric. This restricts the rank of the solution. is set to min(dim(\code{y})-1) by default in "v1".
#' @param lambda Numeric. Nuclear-norm regularization parameter. Controls the low-rank property of the solution
#' to the matrix completion problem. By default, it is determined at the scaling step. If set to zero
#' the algorithm reverts to "hardImputation", where the convergence will be slower. Applicable to "v1" only.
#' @param thresh Numeric. Convergence threshold. Set to 1e-05, by default. Applicable to "v1" only.
#' @param maxit Numeric. Maximum number of iterations of the algorithm before the algorithm is converged. 100 by default.
#' Applicable to "v1" only.
#' @param trace.it Logical. Prints traces of progress of the algorithm.
#' Applicable to "v1" only.
#' @param warm.start List. A SVD object can be used to initialize the algorithm instead of random initialization.
#' Applicable to "v1" only.
#' @param final.svd  Logical. Shall final SVD object be saved?
#' The solutions to the matrix completion problems are computed from U, D and V components of final SVD.
#' Applicable to "v1" only.
#' @param biScale_maxit number of iteration for the scaling algorithm to converge . See \code{scaleData}. You may need to change this
#' parameter only if you're running \code{method=v1}. Applicable to "v1" only.
#' @param gauss_width numeric. The width parameter of the Gaussian distribution to impute the MNAR peptides (features). This the width parameter in the down-shift imputation method.
#' @param gauss_shift numeric. The shift parameter of the Gaussian distribution to impute the MNAR peptides (features). This the width parameter in the down-shift imputation method.
#' @return Missing values are imputed by low-rank approximation of the input matrix. If input is a numeric matrix,
#' a numeric matrix of identical dimensions is returned.
#'
#'
#' @examples
#' data(pxd010943)
#' y <- log2(data.matrix(pxd010943))
#' group <- gsub("_[1234]","", colnames(y))
#' yimp <- msImpute(y, method="v2-mnar", group=group, max.rank=2)
#' @seealso selectFeatures
#' @author Soroor Hediyeh-zadeh
#' @references
#' Hastie, T., Mazumder, R., Lee, J. D., & Zadeh, R. (2015). Matrix completion and low-rank SVD via fast alternating least squares. The Journal of Machine Learning Research, 16(1), 3367-3402.
#' @references
#' Hediyeh-zadeh, S., Webb, A. I., & Davis, M. J. (2020). MSImpute: Imputation of label-free mass spectrometry peptides by low-rank approximation. bioRxiv.
#' @importFrom methods is
#' @export
msImpute <- function(y, method=c("v2-mnar", "v2", "v1"),
                     group = NULL,
                     alpha = 0.2,
		                 relax_min_obs=FALSE,
                     rank.max = NULL, lambda = NULL, thresh = 1e-05,
                     maxit = 100, trace.it = FALSE, warm.start = NULL,
                     final.svd = TRUE, biScale_maxit=20, gauss_width = 0.3, gauss_shift = 1.8) {

  method <- match.arg(method, c("v2-mnar","v2", "v1"))


  if(any(is.nan(y) | is.infinite(y))) stop("Inf or NaN values encountered.")
  
  if(relax_min_obs | any(rowSums(!is.na(y)) <= 3)) {
	  
	  stop("Peptides with excessive NAs are detected. Please revisit your fitering step (at least 4 non-missing measurements are required for any peptide) or set relax_min_obs=TRUE.")
  }
  else if(!relax_min_obs & any(rowSums(!is.na(y)) <= 3)){
	  critical_obs <- which(rowSums(!is.na(y)) <= 3)
	  message("Features with less than 4 non-missing measurements detected. These will be treated as MNAR.")
  }else{
    critical_obs <- NULL
  }
  
  if(any(y < 0, na.rm = TRUE)){
    warning("Negative values encountered in imputed data. Please consider revising filtering and/or normalisation steps.")
  }


  if(!is.null(critical_obs)){
	  y_critical_obs <- y[critical_obs,, drop=FALSE]
    y <- y[-critical_obs,, drop=FALSE]
  }

  if(method=="v1"){
    message(paste("Running msImpute version", method))
    
    yimp <- scaleData(y, maxit = biScale_maxit)
    yimp <- msImputev1(yimp,
                       rank.max = rank.max, lambda = lambda, thresh = thresh,
                       maxit = maxit, trace.it = trace.it, warm.start = warm.start,
                       final.svd = final.svd)
  }else{
    # message(paste("Running msImpute version 2", method))
    message("Running msImpute version 2")
    message("Estimate distribution under MAR assumption")

    rank.max <- ifelse(is.null(rank.max), ceiling(erank(y)) , rank.max)
    yimp <- msImputev1(y, rank.max = rank.max , lambda = estimateLambda(y, rank = rank.max)) #
    if (method == "v2-mnar"){
      message(paste("Compute barycenter of MAR and NMAR distributions", method))
      if (is.null(group)) stop("Please specify the 'group' argument. This is required for the 'v2-mnar' method.")
      ygauss <- gaussimpute(y, width = gauss_width, shift = gauss_shift)
      yimp <- l2bary(y=y, ygauss = ygauss, yerank = yimp, group = group, a=alpha)

    }



  }

  yimp[!is.na(y)] <- y[!is.na(y)]
  if (!is.null(critical_obs)){
	  yimp_critical_obs <- gaussimpute(y_critical_obs, width = gauss_width, shift = gauss_shift)
	  yimp_critical_obs[!is.na(y_critical_obs)] <- y_critical_obs[!is.na(y_critical_obs)]
	  yimp <- rbind(yimp,yimp_critical_obs)
  }


  
  return(yimp)


}


#' @importFrom methods is
#' @keywords internal
msImputev1 <- function(object, rank.max = NULL, lambda = NULL, thresh = 1e-05,
                     maxit = 100, trace.it = FALSE, warm.start = NULL, final.svd = TRUE) {
  # data scaled by biScale
  if(is(object,"list")) {
    x <- object$E
    xnas <- object$E.scaled
  }

  # data is not scaled by biscale
  if(is(object, "matrix")) {
    xnas <- x <- object
    #warning("Input is not scaled. Data scaling is recommended for msImpute optimal performance.")
  }
  # MAList object
  # or \code{MAList} object from \link{limma}
  # if(is(object,"MAList")) x <- object$E

  if(any(is.nan(x) | is.infinite(x))) stop("Inf or NaN values encountered.")
  #if(any(rowSums(!is.na(x)) <= 3)) stop("Peptides with excessive NAs are detected. Please revisit your fitering step. At least 4 non-missing measurements are required for any peptide.")
  if(any(x < 0, na.rm = TRUE)){
    warning("Negative values encountered in imputed data. Please consider revising filtering and/or normalisation steps.")
  }
  if(is.null(rank.max)) rank.max <- min(dim(x) - 1)
  message(paste("rank is", rank.max))
  message("computing lambda0 ...")
  if(is.null(lambda)) lambda <- softImpute::lambda0(xnas)
  message(paste("lambda0 is", lambda))
  message("fit the low-rank model ...")
  fit <- softImpute::softImpute(x, rank.max=rank.max, lambda=lambda,
                                type = "als", thresh = thresh,
                                maxit = maxit, trace.it = trace.it,
                                warm.start = warm.start, final.svd = final.svd)
  message("model fitted. \nImputting missing entries ...")
  ximp <- softImpute::complete(x, fit)
  message("Imputation completed")  # need to define a print method for final rank model fitted

  return(ximp)
  #
  #   if(is(object,"MAList")) {
  #     object$E <- ximp
  #     return(object)
  #   }else{
  #     return(ximp)
  #   }


}

#' @keywords internal
eigenpdf <- function(y, rank=NULL){
  s <- softImpute::softImpute(y, rank.max = ifelse(!is.null(rank), rank, min(dim(y)-1)), lambda =0)$d
  return(s/sum(abs(s)))
}


#' @importFrom stats var sd
#' @keywords internal
estimateS0 <- function(y, rank=NULL){
  set.seed(123)
  s0 <- vector(length = 100L)
  for(i in seq_len(100)){
    s0[i] <- var(eigenpdf(y, rank=rank))
  }
  return(list("s0" = mean(s0), "s0.1sd"= (mean(s0) + sd(s0))))
}

#' @keywords internal
erank <- function(y) {
  P <- eigenpdf(y, rank = NULL)
  return(exp(-sum(P*log(P)))) # shannon entropy
}


#' @keywords internal
estimateLambda <- function(y, rank=NULL) mean(matrixStats::colSds(y, na.rm = TRUE))/estimateS0(y, rank=rank)$"s0.1sd"


#' @importFrom stats quantile
#' @keywords internal
l2bary <- function(y, ygauss, yerank, group, a=0.2){

  pepSds <- matrixStats::rowSds(y, na.rm = TRUE)
  pepMeans <- rowMeans(y, na.rm = TRUE)
  pepCVs <- pepSds/pepMeans
  CV_cutoff <- min(0.2, median(pepCVs))
  varq75 <- quantile(pepSds, p = 0.75, na.rm=TRUE)
  #varq75 <- mean(pepVars)
  EBM <- ebm(y, group)

  # if entropy is nan and variance is low, it is most likely detection limit missing
  # w1 <- ifelse(is.nan(EBM) & (pepCVs < CV_cutoff), 1-a, a)
  w1 <- ifelse(is.nan(EBM), 1-a, a)
  w2 <- 1-w1

  yl2 <- list()
  for(j in colnames(y)){
    yl2[[j]] <- rowSums(cbind(w1*ygauss[,j], w2*yerank[,j]))
  }

  yl2 <- do.call(cbind, yl2)
  yl2[!is.na(y)] <- y[!is.na(y)]
  return(yl2)


}

#' @keywords internal
gaussimpute <- function(x, width=0.3, shift=1.8) {
  # distributions are induced by measured values in each sample
  data.mean <- colMeans(x, na.rm = TRUE)
  data.sd <- matrixStats::colSds(x, na.rm = TRUE)
  n <- nrow(x)
  z <- mvtnorm::rmvnorm(n, mean = data.mean - shift*data.sd , sigma = diag(data.sd*width))
  x[is.na(x)] <- z[is.na(x)]
  return(x)
}


