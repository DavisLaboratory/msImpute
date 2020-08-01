#' Metrics for the assessment of post-imputation structural preservation
#'
#' For an imputed dataset, it computes within phenotype/experimental condition similarity (i.e. preservation of local structures),
#' between phenotype distances (preservation of global structures), and the Gromov-Wasserstein (GW) distance between original and
#' imputed data.
#'
#' @param x numeric matrix. An imputed data matrix.
#' @param group factor. A vector of biological groups, experimental conditions or phenotypes (e.g. control, treatment).
#' @param xna numeric matrix. Data matrix with missing values (i.e. the original intensity matrix with NAs)
#'
#' @details For each group of experimental conditions (e.g. treatment and control), the group centroid is calculated as the average
#' of observed peptide intensities. Withinness for each group is computed as sum of the squared distances between samples in that group and
#' the group centroid. Betweenness is computed as sum of the squared distances between group centroids.
#' When comparing imputation approaches, the optimal imputation strategy should minimize the within
#' group distances, hence smaller withinness, and maximizes between group distances, hence larger betweenness.
#' The GW metric considers preservation of both local and global structures simultaneously. A small GW distance suggests that
#' imputation has introduced small distortions to global and local structures overall, whereas a large distance implies significant
#' distortions. When comparing two or more imputation methods, the optimal method is the method with smallest GW distance.
#' To compute the GW distance, the missing values in each column of \code{xna} are replaced by mean of observed values in that column.
#' This is equivalent to imputation by KNN, where k is set to the total number of identified peptides (i.e. number of rows in the input matrix).
#' GW distance estimation requires \code{python}. See example.
#' All metrics are on log scale.
#'
#'
#' @return list of three metrics: withinness (sum of squared distances within a phenotype group),
#' betweenness (sum of squared distances between the phenotypes), and gromov-wasserstein distance (if \code{xna} is not NULL).
#' All metrics are on log scale.
#'
#'
#' @examples
#' # To compute the GW distance you need to have python installed
#' # then install the reticulate R package from CRAN
#' # install.packages("reticulate")
#' library(reticulate)
#' # create a virtual environment
#' virtualenv_create('r-reticulate')
#' py_available() # if this returns TRUE, you've access to python from R.
#' # See reticulate if you need to troubleshoot
#' # install scipy and POT python modules in this virtual environment
#' virtualenv_install("msImpute-reticulate","scipy")
#' # if this runs successfully, the installation has been successful:
#' scipy <- import("scipy")
#' # You can now run the computeStructuralMetrics() function to compute GW distance.
#' # This setup should only be done for the first use. For all subsequent usages
#' # load the virtual environment that you've created using:
#' library(reticulate)
#' use_virtualenv("msImpute-reticulate")
#' # you can then run the computeStructuralMetrics() function.
#' # Note that the reticulate package should be loaded before loading msImpute.
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
#' y <- xna
#' xna <- scaleData(xna)
#' xcomplete <- msImpute(object=xna)
#' G <- as.factor(sample(1:5, 100, replace = TRUE))
#' computeStructuralMetrics(xcomplete, G, y)
#' @export
computeStructuralMetrics <- function(x, group, xna = NULL){
  out <- list(withinness = log(withinness(x, group)),
       betweenness = log(betweenness(x,group)))

  if(!is.null(xna)){
    GW <- gromov_wasserstein(xna, x)
    out[['gw_dist']] <- GW[[2]]$gw_dist
  }
  return(out)
}



#' @export
withinness <- function(x, class_label){
  within_class_dist <- list()
  for(class in class_label){
  centroid <- colMeans(t(x[,class_label==class]))
  within_class_dist[class] <- sum(as.matrix(pdist::pdist(t(x[,class_label==class]), centroid))^2)
  }
  return(unlist(within_class_dist))
}

#' @export
betweenness <- function(x, class_label){
  centroids <- sapply(unique(class_label), FUN=function(class){
    colMeans(t(x[,class_label==class]))
  })

  return(sum(dist(centroids)^2))
  #return(lsa::cosine(centroids))

}


#' @export
gromov_wasserstein <- function(xna, ximputed){
  reticulate::source_python("gw.py")
  xna <- apply(xna, 2, FUN=function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)})
  return(gw(t(xna), t(ximputed), ncol(xna)))
}
