#' Fills missing values by Peptide Identity Propagation (PIP)
#'
#' Peptide identity (sequence and charge) is propagated from MS-MS or PASEF identified features in evidence.txt to
#' MS1 features in allPeptides.txt that are detected but not identified. A confidence score (probability)
#' is assigned to every propagation. The confidence scores can be used as observation-level weights
#' in \code{limma::lmFit} to account for uncertainty in inferred peptide intensity values.
#'
#' @details
#' Data completeness is maximised by Peptide Identity Propagation (PIP) from runs where
#' a peptide is identified by MSMS or PASEF to runs where peptide is not fragmented
#' (hence MS2 information is not available), but is detected at the MS1 level. \code{mspip} reports a
#' confidence score for each peptide that was identified by PIP. The intensity values of PIP peptides
#' can be used to reduce missing values, while the reported confidence scores can be used to
#' weight the contribution of these peptide intensity values to variance estimation in linear models fitted in
#' \code{limma}.
#'
#' @param path_txt character. The path to MaxQuant \code{txt} directory
#' @param k numeric. The \code{k} nearest neighbors to be used for identity propagation. default to 10.
#' @param thresh numeric. The uncertainty threshold for calling a Identity Transfer as confident. Sequence to peptide
#' feature assignments with confidence score (probability) above a threshold (specified by \code{thresh}) are
#' considered as confident assignments.The rest of the assignments are discarded and not reported in the output.
#' @param skip_weights logical. If TRUE, the propagation confidence scores are also reported.
#' The confidence scores can be used as observation-level weights in \code{limma} linear models
#' to improve differential expression testing. default to FALSE.
#' @param tims_ms logical. Is data acquired by TIMS-MS? default to FALSE.
#'
#' @author Soroor Hediyeh-zadeh
#' @seealso evidenceToMatrix
#'
#' @importFrom dplyr anti_join semi_join
#' @importFrom FNN get.knnx
#' @importFrom utils read.delim
#' @export
mspip <- function(path_txt, k = 10, thresh = 0, skip_weights = TRUE, tims_ms = FALSE){

  evidence_path <- list.files(path=path_txt, pattern = "evidence.txt", full.names = TRUE)
  allPeptides_path <- list.files(path=path_txt, pattern = "allPeptides.txt", full.names = TRUE)

  if(!all(file.exists(allPeptides_path, allPeptides_path))) stop("Required MaxQuant tables are not found in the specified directory")

  message("Reading evidence table")
  evidence <- read.delim(evidence_path,
                         header = TRUE,
                         stringsAsFactors = FALSE)

  # create peptide id
  evidence$PeptideID <- paste0(evidence$Modified.sequence, evidence$Charge)

  # keep only the most intense feature?


  message("Reading allPeptides table")
  allPeptides <- read.delim(allPeptides_path,
                           header = TRUE,
                           stringsAsFactors = FALSE)


  message("Extracting unidentified MS1 peptide features")


  ms1_anchors_pasef <- c("Raw.file","Charge", "Intensity",
                         #"Number.of.isotopic.peaks",
                         "Ion.mobility.index")

  ## MSMS types are problematic here. They aren't proper idents though, so all good.
  ms1_anchors_msms <- c("Raw.file","Charge", "Intensity",
                        # "Number.of.isotopic.peaks",
                        "Number.of.scans")

  ms1_anchors <- ms1_anchors_msms
  if(tims_ms) ms1_anchors <- ms1_anchors_pasef



  identified_peptides <- dplyr::semi_join(evidence, allPeptides,
                                          by = ms1_anchors)
  identified_peptides$Raw.file.id <- as.numeric(as.factor(identified_peptides$Raw.file))
  pep_ids <- as.numeric(as.factor(identified_peptides$PeptideID))
  # pep_f <- as.factor(identified_peptides$PeptideID)

  # do we need RT for matching here? not in PASEF
  unidentified_peptides <- dplyr::anti_join(allPeptides, identified_peptides,
                                            by = ms1_anchors)




  attr_msms <- c("Retention.time","Charge","m.z",
                 "Mass",
                 #"Number.of.isotopic.peaks",
                 "Intensity")

  attr_pasef <- c("Retention.time","Charge","m.z",
                 "Mass",
                 #"Number.of.isotopic.peaks",
                 "Ion.mobility.index",
                 "Intensity")
  attrs <- attr_msms
  if(tims_ms) attrs <- attr_pasef

  query_data <- unidentified_peptides
  query_data$Raw.file.id <- as.numeric(as.factor(query_data$Raw.file))

  # message("Computing one-hot encoding of identifications")
  # one_hot_idents_encoding <- model.matrix(~ 0 + pep_f)
  # C1 <- dplyr::bind_cols(identified_peptides[ , # no keep_idents for rows as what to retain idents in same run as query run
  #                                  c("Retention.time",
  #                                    # "Charge",
  #                                    #"m.z",
  #                                    #"Mass",
  #                                    "Raw.file.id")],
  #             as.data.frame(one_hot_idents_encoding))



  transfered_idents <- list()

  message(paste("Propagating Peptide Identities within", k, "nearest neighbors per run"))
  for (run_id in unique(evidence$Raw.file)){
    message(run_id)
    missing_idents <- setdiff(identified_peptides$PeptideID[!identified_peptides$Raw.file %in% run_id & !is.na(identified_peptides$Intensity)],
                              identified_peptides$PeptideID[identified_peptides$Raw.file %in% run_id & !is.na(identified_peptides$Intensity)])

    # run_idents <- unique(identified_peptides$PeptideID[identified_peptides$Raw.file %in% run_id & !is.na(identified_peptides$Intensity)])

    message("Number of missing idents")
    message(length(missing_idents))




    keep_idents <- (identified_peptides$PeptideID %in% missing_idents) & (!identified_peptides$Raw.file %in% run_id)


    # compute width of Random Walk

    # sigma <- matrixStats::rowMedians(FNN::get.knn(identified_peptides[keep_idents,
    #                                                      c("Retention.time","Charge",
    #                                                       "m.z","Mass",
    #                                                       "Mod..peptide.ID",
    #                                                       "Number.of.isotopic.peaks","Intensity")],
    #                            k = 5)$nn.dist)


    # message("sigma")
    # message(sqrt(sigma))





    # C2 <- query_data[query_data$Raw.file %in% run_id, c("Raw.file.id","Retention.time")]
    # one_hot_encoding_query <- matrix(0,nrow(C2), max(pep_ids))
    # C2 <- cbind(C2, one_hot_encoding_query)
    # elutions <- rbind(C1,C2)
    # coelutions <- dbscan::sNN(elutions, k = 5, kt = 5)


    # message("Building sNN graphs")
    # elutions <- identified_peptides[keep_idents , c("Retention.time", "Raw.file.id")]
    # snn_elutions_donor_runs <- dbscan::sNN(elutions, k = 5, kt = 3)
    #
    #
    #
    # coelute_idents <- matrix(pep_ids[keep_idents][snn_elutions_donor_runs$id],
    #                          byrow=FALSE,
    #                          nrow = nrow(snn_elutions_donor_runs$id),
    #                          ncol = ncol(snn_elutions_donor_runs$id))
    #
    # coelute_idents[is.na(coelute_idents)] <- 0

    # coelute_mz <- matrix(identified_peptides$m.z[keep_idents][coelutions$nn.index],
    #                          byrow=FALSE,
    #                          nrow = nrow(coelutions$nn.index),
    #                          ncol = ncol(coelutions$nn.index))
    #
    # coelute_rt <- matrix(identified_peptides$Retention.time[keep_idents][coelutions$nn.index],
    #                      byrow=FALSE,
    #                      nrow = nrow(coelutions$nn.index),
    #                      ncol = ncol(coelutions$nn.index))








    # identifications
    run_prototypes <- identified_peptides[keep_idents, attrs]
    #run_prototypes <- cbind(run_prototypes, coelute_idents)

    ident_labels <- identified_peptides[keep_idents, "PeptideID"]
    prototype_charges <- as.numeric(run_prototypes$Charge)


    # detected features
    query_embedding <- query_data[query_data$Raw.file %in% run_id, attrs]
    query_charge <- as.numeric(query_embedding$Charge)


    ### add coelution for query LC-MS features
    # C1 <- identified_peptides[identified_peptides$Raw.file %in% run_id, c("Raw.file.id","Retention.time")]
    # query_coelutions <- FNN::get.knnx(query_elutions,
    #                                   query_data[query_data$Raw.file %in% run_id, c("Raw.file.id","Retention.time")],
    #                                   k = 5)


    # C2 <- query_data[query_data$Raw.file %in% run_id, c("Raw.file.id","Retention.time")]
    #
    # query_elutions <- rbind(C1, C2)
    # snn_elutions_query <- dbscan::sNN(query_elutions, k = 5, kt = 3)
    #
    # snn_elutions_query_ids <- snn_elutions_query$id[(nrow(C1) + 1):nrow(query_elutions),]
    #
    # # NA indicies or those larger than nrow C1 are unidentified sNN and should be removed
    # snn_elutions_query_ids[snn_elutions_query_ids > nrow(C1) | is.na(snn_elutions_query_ids)] <- NA
    #
    # query_coelute_idents <- matrix(pep_ids[identified_peptides$Raw.file %in% run_id][snn_elutions_query_ids],
    #                          byrow=FALSE,
    #                          nrow = nrow(snn_elutions_query_ids),
    #                          ncol = ncol(snn_elutions_query_ids))
    #
    # query_coelute_idents[is.na(query_coelute_idents)] <- 0

    # query_coelute_mz <- matrix(identified_peptides$m.z[identified_peptides$Raw.file %in% run_id][query_coelutions$nn.index],
    #                          byrow=FALSE,
    #                          nrow = nrow(query_coelutions$nn.index),
    #                          ncol = ncol(query_coelutions$nn.index))
    #
    # query_coelute_rt <- matrix(identified_peptides$Retention.time[identified_peptides$Raw.file %in% run_id][query_coelutions$nn.index],
    #                          byrow=FALSE,
    #                          nrow = nrow(query_coelutions$nn.index),
    #                          ncol = ncol(query_coelutions$nn.index))

    #query_embedding <- cbind(query_embedding, query_coelute_idents)


    message("Number of detected features available for PIP in the run")
    message(nrow(query_embedding))
    # knn_prototypes <- FNN::get.knnx(run_prototypes[, grep("Intensity", colnames(run_prototypes), invert = TRUE)],
    #                                 query_embedding[, grep("Intensity", colnames(query_embedding), invert = TRUE)], k = 10) # nsamples - 1


    message("Propagating identifications")
    knn_prototypes <- FNN::get.knnx(
      query_embedding[, grep("Intensity", colnames(query_embedding), invert = TRUE)],
      run_prototypes[, grep("Intensity", colnames(run_prototypes), invert = TRUE)],
      k = k) # nsamples - 1


    # probs <- exp(-0.5*((knn_prototypes$nn.dist^2))) # i.e. sigma = 1
    # probs <- exp(-0.5*((knn_prototypes$nn.dist^2)/matrixStats::rowSds(knn_prototypes$nn.dist)))

    probs <- exp(-0.5*((knn_prototypes$nn.dist^2)/matrixStats::rowMedians(knn_prototypes$nn.dist)))

    # probs <- exp(-0.5*((knn_prototypes$nn.dist^2)/sigma))



    # ww <- matrix(prototype_charges[knn_prototypes$nn.index], nrow = nrow(probs), ncol = ncol(probs))
    # charge <- matrix(query_charge, nrow = nrow(ww), ncol = ncol(ww), byrow = FALSE)


    ww <- matrix(query_charge[knn_prototypes$nn.index], nrow = nrow(probs), ncol = ncol(probs))
    charge <- matrix(prototype_charges, nrow = nrow(ww), ncol = ncol(ww), byrow = FALSE)

    w <- ifelse(ww==charge, 1, 0)

    wprobs <- w*probs

    p1 <- wprobs
    p2 <- wprobs/rowSums(probs)
    p3 <- wprobs/rowSums(wprobs)


    normalised_probs <- p3

    if(sum(!complete.cases(normalised_probs)) > 0 ) {
      message("Warning: No MS1 feature was found for some identifications.You may wish to increase k.")
    }

    valid_features <- rowSums(is.finite(normalised_probs)) > 1

    normalised_probs <- normalised_probs[valid_features,]
    nn_indices <- knn_prototypes$nn.index[valid_features,]

    idxs <- apply(normalised_probs, 1, FUN= function(x) {
      z <- logical(length(x)); z[which.max(x)] <- TRUE; return(z)
    })

    idxs <- matrix(as.vector(idxs), nrow = nrow(normalised_probs),
                   ncol = ncol(normalised_probs),
                   byrow = FALSE)
    max_probs <- t(normalised_probs)[idxs]

    query_max_probs <- t(nn_indices)[idxs]

    df_query_idents <- cbind(
      Raw.file = run_id,
      query_embedding[query_max_probs, grep("[1-9]", colnames(query_embedding), invert = TRUE)],
      data.frame(probability = max_probs, PeptideID = ident_labels[valid_features])
    )

    rownames(df_query_idents) <- NULL
    transfered_idents[[run_id]] <- df_query_idents

  }

  transfered_idents <- do.call(rbind, transfered_idents)
  message(paste("Discarding", sum(!(transfered_idents$probability > thresh)),
                "low-confidence PIPs at threshold", thresh))
  if(skip_weights){
    evidence_pip <- rbind(evidence[,c("Raw.file","PeptideID", "Intensity")],
                          transfered_idents[transfered_idents$probability > thresh,
                                            c("Raw.file","PeptideID", "Intensity")])
  }else{
    evidence_pip <- rbind(
      cbind(evidence[,c("Raw.file","PeptideID", "Intensity")], weight = 1),
      cbind(transfered_idents[transfered_idents$probability > thresh,
                              c("Raw.file","PeptideID", "Intensity")],
            weight = transfered_idents$probability[transfered_idents$probability > thresh])
                          )


    meta_attrs <- c( "PeptideID", "Sequence", "Length", "Modifications",
                     "Modified.sequence",
                     "Leading.razor.protein","Gene.Names", "Protein.Names",
                     "Charge")
    evidence_colnames <- tolower(colnames(evidence))

    genes <- evidence[,match(tolower(meta_attrs), evidence_colnames)]
    genes <- genes[!duplicated(genes),]
    evidence_pip <- cbind(evidence_pip, genes[match(evidence_pip$PeptideID, genes$PeptideID), grep("PeptideID", colnames(genes), invert=TRUE)])
  }
  message("PIP completed")
  return(evidence_pip)

}
