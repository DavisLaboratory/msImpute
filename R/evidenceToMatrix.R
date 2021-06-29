#' Creates intensity matrix from tabular data in evidence table of MaxQuant
#'
#' Every \code{Modified sequence} - \code{Charge} is considered as a precursor feature.
#' Only the feature with maximum intensity is retained. The columns are run names, the rows
#' are peptide ids (in the \code{Modified.sequence_Charge} format)
#'
#' @param evidence data.frame. The evidence table read from evidence.txt, or data.frame created by \code{mspip}.
#' @param run_id character. The name of the column of evidence containing the run/raw file name.
#' These form the columns of the intensity data matrix.
#' @param peptide_id character. The name of the column of evidence containing the peptide ids.
#' These form the rows of the intensity data matrix.
#' @param return_EList logical. If TRUE, returns a \code{EListRaw} object storing both the
#' intensity data matrix and observation-level weights from
#' \code{mspip} (propagation confidence score), otherwise returns a matrix.
#' @param weights character. The name of the column of evidence containing weights from \code{mspip}. default to NULL.
#' Set this to "weight" if you want the weights from PIP stored in the \code{weights} slot of the \code{EListRaw} object.
#'
#'
#' @return a numeric matrix of intensity data, or a \code{EListRaw} object containing
#' such data and observation-level weights from \code{mspip}.
#'
#' @details The \code{EListRaw} object created by the function is intended to bridge \code{msImpute} and statistical
#' methods of \code{limma}. The object can be passed to \code{normalizeBetweenArrays} for normalisation, which can then
#' be passed to \code{lmFit} and \code{eBayes} for fitting linear models per peptide and Empirical Bayes moderation of t-statistics
#' respectively. The \code{weights} slot is recognized by \code{lmFit}, which incorporates the uncertainty in intensity values
#' inferred by PIP into the test statistic.
#' The function is also a generic tool to create a matrix or \code{limma}-compatible objects from the evidence table of MaxQuant.
#'
#' @importFrom stats aggregate
#' @importFrom tidyr spread
#' @importFrom stats na.pass complete.cases
#' @importFrom methods new
#' @seealso mspip
#' @export
#' @author Soroor Hediyeh-zadeh
evidenceToMatrix <- function(evidence, run_id = "Raw.file", peptide_id = "PeptideID",
                             return_EList = FALSE, weights = NULL){



  y <- aggregate(evidence[,"Intensity"] ~ evidence[, run_id] + evidence[, peptide_id],
                 FUN = function(x) max(x, na.rm=TRUE),
                 na.action = na.pass)

  colnames(y) <- c(run_id, peptide_id, "Intensity")
  y[y==-Inf] <- NA

  E <- tidyr::spread(y, key = 1, value = 3)

  rownames(E) <- E[,1]
  E <- E[,-1]
  #E[E == -Inf] <- NA

  E <- data.matrix(E)

  if(return_EList){

    meta_attrs <- c( peptide_id, "Sequence", "Length", "Modifications",
                     "Modified.sequence",
                     "Leading.razor.protein","Gene.Names", "Protein.Names",
                     "Charge")
    evidence_colnames <- tolower(colnames(evidence))

    # genes <- evidence[,match(tolower(meta_attrs), evidence_colnames)]
    genes <- evidence[, evidence_colnames %in% tolower(meta_attrs)]
    genes <- genes[!duplicated(genes),]
    genes <- genes[match(rownames(E), genes[,peptide_id]),]


    if(!is.null(weights)){
      if (!weights %in% colnames(evidence)) {
        message("No weight column in the input. Returning an EList without the weights slot")
        return(new("EListRaw", list(E=E, genes = genes)))
      } else{
        idx <- match(paste0(y[,run_id], y[,peptide_id],y[,"Intensity"]),
                     paste0(evidence[,run_id], evidence[,peptide_id], evidence[,"Intensity"])
        )
        w <- evidence[idx, c(run_id, peptide_id, "weight")]
        weights <- tidyr::spread(w, key = 1, value = 3)
        rownames(weights) <- weights[,1]
        weights <- weights[,-1]
        weights[is.na(weights)] <- 0 # when pip idents are filtered, NAs will appear in weight matrix.

        return(new("EListRaw", list(E=E, weights=weights, genes = genes)))
      }

    } else{
      return(new("EListRaw", list(E=E, genes = genes)))
    }
  } else {
    return(E)
  }

}
