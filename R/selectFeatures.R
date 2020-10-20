#' Select features for MAR/MNAR pattern examination
#'
#' Two methods are provided to identify features (peptides or proteins) that can be informative of missing patterns.
#' Method \code{hvp} fits a linear model to peptide dropout rate (proportion of samples were peptide is missing)
#' against peptide abundance (average log2-intensity). Method \code{emb} is a information theoretic approach to
#' identify missing patterns. It quantifies the heterogeneity (entropy) of missing patterns per
#' biological (experimental group). This is the default method.
#'
#' @details
#' In general, the presence of group-wise (structured) blocks of missing values,
#' where peptides are missing in one experimental group can indicate MNAR, whereas if
#' such patterns are absent (or missingness is uniform across the samples), peptides are likely MAR.
#' In the presence of MNAR, left-censored MNAR imputation methods should
#' be chosen. Two methods are provided to explore missing patterns: \code{method=hvp} identifies top \code{n_features}
#' peptides with high average expression that also have high dropout rate, defined as the proportion of samples where
#' peptide is missing. Peptides with high (potentially) biological dropouts are marked in the \code{hvp} column in the
#' output dataframe. This method does not use any information about experimental conditions (i.e. group).
#' Another approach to explore and quantify missing patterns is by looking at how homogeneous or heterogeneous
#' missing patterns are in each experimental group. This is done by computing entropy of distribution of observed values.
#' This is the default and recommended method for \code{selectFeatures}. Entropy is reported in \code{EBM} column
#' of the output. A \code{NaN} EBM indicates peptide is missing at least in one experimental group. Features set to
#' \code{TRUE} in \code{msImpute_feature} column are the features selected by the selected method. Users are encouraged
#' to use the EBM metric to find informative features, hence why the \code{group} argument is required.
#'
#'
#'
#' @param x Numeric matrix giving log-intensity where missing values are denoted by NA.
#' Rows are peptides, columns are samples.
#' @param method character. What method should be used to find features? options include \code{method='hvp'} and \code{method='ebm'}
#' @param group character or factor vector specifying biological (experimental) group e.g. control, treatment, WT, KO
#' @param n_features Numeric, number of features with high dropout rate. 500 by default. Applicable if \code{method="hvp"}.
#' @param suppress_plot Logical show plot of dropouts vs abundances. Default to TRUE. Applicable if \code{method="hvp"}.
#'
#' @return A data frame with a logical column denoting the selected features
#'
#' @examples
#' data(pxd007959)
#' group <- pxd007959$samples$group
#' y <- data.matrix(pxd007959$y)
#' y <- log2(y)
#' hdp <- selectFeatures(y, method="ebm", group = group)
#' # construct matrix M to capture missing entries
#' M <- ifelse(is.na(y),1,0)
#' M <- M[hdp$msImpute_feature,]
#' # plot a heatmap of missingness patterns for the selected peptides
#' require(ComplexHeatmap)
#' hm <- Heatmap(M,
#' column_title = "dropout pattern, columns ordered by dropout similarity",
#'               name = "Intensity",
#'               col = c("#8FBC8F", "#FFEFDB"),
#'               show_row_names = FALSE,
#'               show_column_names = TRUE,
#'               cluster_rows = TRUE,
#'               cluster_columns = TRUE,
#'               show_column_dend = TRUE,
#'               show_row_dend = FALSE,
#'               row_names_gp =  gpar(fontsize = 7),
#'               column_names_gp = gpar(fontsize = 8),
#'               heatmap_legend_param = list(#direction = "horizontal",
#'               heatmap_legend_side = "bottom",
#'               labels = c("observed","missing"),
#'               legend_width = unit(6, "cm")),
#'          )
#' hm <- draw(hm, heatmap_legend_side = "left")
#' @author Soroor Hediyeh-zadeh
#' @seealso msImpute
#' @references
#' Hediyeh-zadeh, S., Webb, A. I., & Davis, M. J. (2020). MSImpute: Imputation of label-free mass spectrometry peptides by low-rank approximation. bioRxiv.
#' @importFrom stats lm residuals
#' @importFrom methods is
#' @importFrom graphics abline plot
#' @export
selectFeatures <- function(x, method=c("ebm","hvp"), group, n_features=500, suppress_plot = TRUE) {

  if(is.null(rownames(x))) stop("No row names in input. Please provide input with named rows.")
  if(any(is.nan(x) | is.infinite(x))) stop("Inf or NaN values encountered.")

  AveExpr <- rowMeans(x, na.rm = TRUE)
  dropout <- rowMeans(is.na(x))

  linear_fit <- lm(dropout ~ AveExpr)
  resids <- residuals(linear_fit)
  lin_res_o <- order(resids, decreasing = TRUE)

  # Entropy of batch mixing----
  EBM <- ebm(x=x,group=group)


  # default method is ebm
  method <- match.arg(method, c("ebm","hvp"))

  if(!suppress_plot & method=="hvp"){
    cols <- rep("#3E71A8", length(resids))
    cols[lin_res_o[seq_len(n_features)]] <- "#DE1A1A"
    plot(x = AveExpr, y = dropout, pch = 16,
         cex = 0.5, col = cols, main = paste("Top ",n_features," high droupout peptides", sep =""))
    abline(linear_fit)
  }

  hdrp <- data.frame(name = rownames(x), AveExpr = AveExpr, dropout = dropout,
                         residual = resids, hvp=FALSE, EBM=EBM, msImpute_feature=FALSE)


  hdrp$hvp[lin_res_o[seq_len(n_features)]] <- TRUE


  if(method=="hvp"){
    hdrp$msImpute_feature[lin_res_o[seq_len(n_features)]] <- TRUE
  }

  if(method=="ebm"){
    if(all(!is.nan(EBM))){
      message("No NaN EBMs detected. Peptides are missing evenly across samples.")
      message("Switchted to 'hvp' method as final msImpute features")
      hdrp$msImpute_feature[lin_res_o[seq_len(n_features)]] <- TRUE
    }else{
      hdrp$msImpute_feature[is.nan(EBM)] <- TRUE
    }

  }

  hdrp <- data.table::as.data.table(hdrp)

  return(hdrp)

}

#' @keywords internal
ebm <- function(x, group){
  M <- ifelse(is.na(x), 1,0)
  P <- list()
  for(i in unique(group)){
    P[[i]] <- rowMeans(M[,group==i]==0)*log(rowMeans(M[,group==i]==0)) # i.e. number observed entries per group
  }

  Pmat <- do.call(cbind, P)
  return(-rowSums(Pmat))
}
