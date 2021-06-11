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
#' @param return_MAList logical. If TRUE, returns a \code{MAList} object storing both the
#' intensity data matrix and observation-level weights from
#' \code{mspip} (propagation confidence score), otherwise returns a matrix.
#'
#'
#' @return a numeric matrix of intensity data, or a \code{MAList} object containing
#' such data and observation-level weights from \code{mspip}
#'
#' @importFrom stats aggregate
#' @importFrom tidyr spread
#' @importFrom stats na.pass
#' @seealso mspip
#' @export
#' @author Soroor Hediyeh-zadeh
evidenceToMatrix <- function(evidence, run_id = "Raw.file", peptide_id = "PeptideID",
                             return_MAList = FALSE){



  y <- aggregate(evidence[,"Intensity"] ~ evidence[, run_id] + evidence[, peptide_id],
                 FUN = function(x) max(x, na.rm=TRUE),
                 na.action = na.pass)

  y <- tidyr::spread(y, key = 1, value = 3)
  rownames(y) <- y[,1]
  y <- y[,-1]
  y[y == -Inf] <- NA

  y <- data.matrix(y)

}
