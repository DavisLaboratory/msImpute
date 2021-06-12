#' Creates intensity matrix from tabular data in evidence table of MaxQuant
#'
#' Every \code{Modified sequence} - \code{Charge} is considered as a precursor feature.
#' Only the feature with maximum intensity is retained. The columns are run names, the rows
#' are peptide ids (in the \code{Modified.sequence_Charge} format)
#'
#' @param evidence data.frame. The evidence table read from evidence.txt, or \code{mspip} output
#' @param run_id character. The name of the column of evidence containing the run/raw file name.
#' These form the columns of the intensity data matrix.
#' @param peptide_id character. The name of the column of evidence containing the peptide ids.
#' These form the rows of the intensity data matrix.
#' @param return_EList logical. If TRUE, returns a \code{EListRaw} object storing both the
#' intensity data matrix and observation-level weights from
#' \code{mspip} (propagation confidence score), otherwise returns a matrix.
#'
#'
#' @return a numeric matrix of intensity data, or a \code{EListRaw} object containing
#' such data and observation-level weights from \code{mspip}
#'
#' @importFrom stats aggregate
#' @importFrom tidyr spread
#' @importFrom stats na.pass complete.cases
#' @importFrom methods new
#' @seealso mspip
#' @export
#' @author Soroor Hediyeh-zadeh
evidenceToMatrix <- function(evidence, run_id = "Raw.file", peptide_id = "PeptideID",
                             return_EList = FALSE){



  y <- aggregate(evidence[,"Intensity"] ~ evidence[, run_id] + evidence[, peptide_id],
                 FUN = function(x) max(x, na.rm=TRUE),
                 na.action = na.pass)

  colnames(y) <- c(run_id, peptide_id, "Intensity")


  E <- tidyr::spread(y, key = 1, value = 3)

  rownames(E) <- E[,1]
  E <- E[,-1]
  E[E == -Inf] <- NA

  E <- data.matrix(E)

  if(return_EList){
    if (!"weight" %in% colnames(evidence)) stop("No weight column in the input.")
    idx <- match(paste0(y[,run_id], y[,peptide_id],y[,"Intensity"]),
                 paste0(evidence[,run_id], evidence[,peptide_id], evidence[,"Intensity"])
                 )
    w <- evidence[idx, c(run_id, peptide_id, "weight")]
    weights <- tidyr::spread(w, key = 1, value = 3)
    rownames(weights) <- weights[,1]
    weights <- weights[,-1]
    weights[is.na(weights)] <- 0 # when pip idents are filtered, NAs will appear in weight matrix.

    return(new("EList", list(E=E, weights=weights)))
  } else {
    return(E)
  }

}