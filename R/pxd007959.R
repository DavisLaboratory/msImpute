#' Processed peptide intensity matrix and experimental design table from PXD007959 study
#'
#' Extracellular vesicles isolated from the descending colon of pediatric patients with inflammatory bowel disease
#' and control patients. Characterizes the proteomic profile of extracellular vesicles isolated from the descending colon
#' of pediatric patients with inflammatory bowel disease and control participants. This object contains data from peptide.txt
#' table output by MaxQuant. Rows are Modified Peptide IDs. Charge state variations are treated as distinct peptide species.
#' Reverse complements and contaminant peptides are discarded. Peptides with more than 4 observed intensity values are retained.
#' Additionally, qualified peptides are required to map uniquely to proteins.
#' Two of the samples with missing group annotation were excluded.
#' The peptide.txt and experimentalDesignTemplate files can be downloaded as RDS object from \url{https://github.com/soroorh/proteomicscasestudies}.
#' Code for data processing is provided in package vignette.
#'
#' @format A list of two: samples (data frame of sample descriptions), and y (numeric matrix of peptide intensity values)
#' @references
#' Zhang X, Deeke SA, Ning Z, Starr AE, Butcher J, Li J, Mayne J, Cheng K, Liao B, Li L, Singleton R, Mack D, Stintzi A, Figeys D, Metaproteomics reveals associations between microbiome and intestinal extracellular vesicle proteins in pediatric inflammatory bowel disease. Nat Commun, 9(1):2873(2018)
#' @source \url{http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD007959}
"pxd007959"
