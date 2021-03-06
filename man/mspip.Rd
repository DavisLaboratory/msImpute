% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspip.R
\name{mspip}
\alias{mspip}
\title{Fills missing values by Peptide Identity Propagation (PIP)}
\usage{
mspip(
  path_txt,
  k = 10,
  thresh = 0,
  skip_weights = TRUE,
  tims_ms = FALSE,
  group_restriction = NULL,
  nlandmarks = 50
)
}
\arguments{
\item{path_txt}{character. The path to MaxQuant \code{txt} directory}

\item{k}{numeric. The \code{k} nearest neighbors to be used for identity propagation. default to 10.}

\item{thresh}{numeric. The uncertainty threshold for calling a Identity Transfer as confident. Sequence to peptide
feature assignments with confidence score (probability) above a threshold (specified by \code{thresh}) are
considered as confident assignments.The rest of the assignments are discarded and not reported in the output.}

\item{skip_weights}{logical. If TRUE, the propagation confidence scores are also reported.
The confidence scores can be used as observation-level weights in \code{limma} linear models
to improve differential expression testing. default to FALSE.}

\item{tims_ms}{logical. Is data acquired by TIMS-MS? default to FALSE.}

\item{group_restriction}{A data.frame with two columns named Raw.file and group, specifying run file and the (experimental) group to which the run belongs.
Use this option for Unbalanced PIP}

\item{nlandmarks}{numeric. Number of landmark peptides used for measuring neighborhood/coelution similarity. Default to 50.}
}
\description{
Peptide identity (sequence and charge) is propagated from MS-MS or PASEF identified features in evidence.txt to
MS1 features in allPeptides.txt that are detected but not identified. A confidence score (probability)
is assigned to every propagation. The confidence scores can be used as observation-level weights
in \code{limma::lmFit} to account for uncertainty in inferred peptide intensity values.
}
\details{
Data completeness is maximised by Peptide Identity Propagation (PIP) from runs where
a peptide is identified by MSMS or PASEF to runs where peptide is not fragmented
(hence MS2 information is not available), but is detected at the MS1 level. \code{mspip} reports a
confidence score for each peptide that was identified by PIP. The intensity values of PIP peptides
can be used to reduce missing values, while the reported confidence scores can be used to
weight the contribution of these peptide intensity values to variance estimation in linear models fitted in
\code{limma}.
}
\seealso{
evidenceToMatrix
}
\author{
Soroor Hediyeh-zadeh
}
