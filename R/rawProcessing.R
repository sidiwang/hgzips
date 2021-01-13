#' rawProcessing
#'
#' This \code{rawProcessing} function finds the actual (N_ij) and expected counts (E_ij) of each AE - drug pair by implementing the methodology described by DuMouchel (1999); This function outputs N and E in both matrix format and vector format.
#' @import stats
#' @import openEBGM
#'
#' @param rawdata A data frame containing columns named: primaryid, PT (Adverse Events), prod_ai (drug), strat1 (optional), strat2 (optional), ...... stratx
#' @param stratify A logical scalar specifiying whether stratification is to be used (TRUE) or not (FALSE)
#' @param frequency_threshold the minimum frequency of each AE or drug for it to be kept in the dataset
#' @param data.squashing whether to conduct data squashing or not (TRUE or FALSE)
#' @param zeroes A logical scalar specifying if zero counts should be included.
#' @param keep_pts A vector of whole numbers for the number of points to leave unsquashed for each count (N). See the 'openEBGM' details section.
#'
#' @return a list including the following:
#' \itemize{
#' \item \code{N_ij} actual counts of each AE - drug pair in matrix format
#' \item \code{Nij} actual counts of each AE - drug pair in vector format
#' \item \code{E_ij} expected counts of each AE - drug pair in matrix format
#' \item \code{Eij} expected counts of each AE - drug pair in vector format
#' \item \code{processedData} dataset after deleting AE or drug with frequencies less than \code{frequency_threshold}
#' \item \code{drugList} a vector of drug names that are kept in \code{processedData}
#' \item \code{AElist} a vector of adverse event names that are kept in \code{processedData}
#' \item \code{N_ij_squashed} a list of squashed N - E - weight dataset for each drug
#' \item \code{N_ij_list} a list of original N - E - weight (= 1) dataset for each drug
#' \item \code{squashed} squashed \code{processedData} dataset
#' }
#' @details An \code{primaryid} column must be included. Columns must be named as as \code{primaryid, PT, prod_ai, strat1, strat2, ... stratx}. Only variables containing 'strat' (case sensitive) will be used as stratification variables.
#' @export
#' @seealso openEBGM
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency Tables, With an Application to the FDA Spontaneous Reporting System." The American Statistician, 53(3), 177-190.
#' openEBGM R - package
#' @examples
#' primaryid = c(1:10)
#' PT = c("fatigue", "restlessness", "sleepdisorder", "fatigue", "anaemia", "blood creatinine increased", "eosinophilia", "generalised oedema", "hypertension", "nephropathy toxic")
#' prod_ai = c("peginterferon alfa-2a", "peginterferon alfa-2a", "peginterferon alfa-2a", "pantoprazole", "ribavirin", "cyclosporine", "cyclosporine", "cyclosporine", "ribavirin", "flunitrazepam")
#' strat1 = c("M", "M", "M", "F", "F", "M", "M", "M", "F", "F")
#' strat2 = c("Young", "Young", "Young", "Old", "Old", "Old", "Old", "Young", "Young", "Young")
#' dat = data.frame(primaryid = primaryid, PT = PT, prod_ai = prod_ai, strat1 = strat1, strat2 = strat2)
#' result = hgzips::rawProcessing(dat, stratify = TRUE, frequency_threshold = 0, data.squashing = FALSE, zeroes = TRUE, keep_pts = 190000)

############################################################
## raw data processing, generating N_ij, E_ij, Nij, Eij
############################################################

# rawdata has four columns: primaryid, PT, drugname, prod_ai
# input: raw data
# output: two formats of Nij, two formats of Eij


rawProcessing = function(rawdata, stratify = FALSE, frequency_threshold, data.squashing = FALSE, zeroes = FALSE, keep_pts){


  small_tidy = unique(rawdata[, c("primaryid", "PT", "prod_ai")])

  # list of unique PT
  PT_list = unique(small_tidy$PT)
  # list of unique drug
  prod_list = sort(unique(small_tidy$prod_ai))

  # raw N_ij matrix
  N_ij_raw = table(small_tidy$PT, small_tidy$prod_ai)

  # N_ij matrix after removing low frequency drugs and adverse events (threshold = 100)
  rows_to_drop = which(rowSums(N_ij_raw) <= frequency_threshold)
  columns_to_drop = which(colSums(N_ij_raw) <= frequency_threshold)
  if (length(rows_to_drop) != 0){
    if (length(columns_to_drop) == 0){
      N_ij = as.data.frame.matrix(N_ij_raw[-rows_to_drop,])
    } else {
      N_ij = as.data.frame.matrix(N_ij_raw[-rows_to_drop, -columns_to_drop])
    }
  }else{
    if (length(columns_to_drop) == 0){
      N_ij = as.data.frame.matrix(N_ij_raw)
    }else{
      N_ij = as.data.frame.matrix(N_ij_raw[, -columns_to_drop])
    }
  }


  drug_keep = colnames(N_ij)
  AE_keep = rownames(N_ij)

  if (stratify == FALSE){

    # another format
    if (length(rows_to_drop) != 0){
      if (length(columns_to_drop) == 0){
        Nij = as.data.frame(N_ij_raw[-rows_to_drop, ])
      }else{
        Nij = as.data.frame(N_ij_raw[-rows_to_drop, -columns_to_drop])
      }
    }else{
      if (length(columns_to_drop) == 0){
        Nij = as.data.frame(N_ij_raw)
      }else{
        Nij = as.data.frame(N_ij_raw[, -columns_to_drop])
      }
    }
    colnames(Nij) = c("PT", "prod_ai", "frequency")

    # Calculate E_ij
    N_i. = rowSums(N_ij)
    N_.j = colSums(N_ij)
    N_.. = sum(N_ij)

    E_ij = matrix(NA, nrow(N_ij), ncol(N_ij))
    rownames(E_ij) = rownames(N_ij)
    colnames(E_ij) = colnames(N_ij)
    for (j in 1 : ncol(N_ij)){
      E_ij[ , j] = N_i. * N_.j[j] / N_..
    }
  } else {
    small_tidy = rawdata[rownames(unique(rawdata[, c("primaryid", "PT", "prod_ai")])),]
    small_tidy = subset(small_tidy, PT %in% AE_keep & prod_ai %in% drug_keep)
    cat_strat = unique(select(small_tidy, matches("strat")))
    cat_num = nrow(cat_strat)
    cat_col = ncol(cat_strat)
    strat_names = colnames(cat_strat)

    for (i in 1:cat_num){
      assign(paste0("N_ij_", i), small_tidy)
      for (j in 1:cat_col){
        assign(paste0("N_ij_", i), subset(get(paste0("N_ij_", i)), get(paste0("N_ij_", i))[, strat_names[j]] == cat_strat[i, j]))
      }
    }

    for (i in 1:cat_num){

      data = get(paste0("N_ij_", i))
      data1 = table(data$PT, data$prod_ai)
      data2 = matrix(NA, nrow(N_ij), ncol(N_ij))
      colnames(data2) = colnames(N_ij)
      rownames(data2) = rownames(N_ij)
      data2[rownames(data2) %in% rownames(data1), colnames(data2) %in% colnames(data1)] <- data1[rownames(data1) %in% rownames(data2), colnames(data1) %in% colnames(data2)]
      data2[is.na(data2)] = 0
      assign(paste0("N_strat", i), data2)

    }

    N_ij_strat = matrix(0, nrow(N_ij), ncol(N_ij))
    for (i in 1:cat_num){
      N_ij_strat = N_ij_strat + get(paste0("N_strat", i))
    }

    for (i in 1:cat_num){
      assign(paste0("N_i.", i), rowSums(get(paste0("N_strat", i))))
      assign(paste0("N_.j", i), colSums(get(paste0("N_strat", i))))
      assign(paste0("N_..k", i), sum(get(paste0("N_strat", i))))
    }

    E_ij_strat = matrix(0, nrow(N_ij), ncol(N_ij))
    colnames(E_ij_strat) = colnames(N_ij)
    rownames(E_ij_strat) = rownames(N_ij)

    for (i in 1:cat_num){
      for (j in 1 : ncol(N_ij)){
        E_ij_strat[, j] = E_ij_strat[, j] + get(paste0("N_i.", i)) * get(paste0("N_.j", i))[j] / get(paste0("N_..k", i))
      }
    }

    E_ij = E_ij_strat

    N_ij = N_ij_strat

    # another format
    Nij = as.data.frame(as.table(N_ij))
    colnames(Nij) = c("PT", "prod_ai", "frequency")

  }



  Eij = as.data.frame(as.table(E_ij))
  colnames(Eij) = c("PT", "prod_ai", "baseline")

  E_ij = as.data.frame(E_ij)

  if (data.squashing == FALSE){

    N_E_list = list()
    for (j in (1:ncol(N_ij))){
      jj = as.data.frame(cbind(N_ij[, j], E_ij[, j], rep(1, nrow(N_ij))))
      colnames(jj) = c("N", "E", "weight")
      N_E_list[[j]] = jj
    }

    NnE <- list("N_ij" = N_ij, "Nij" = Nij, "E_ij" = E_ij, "Eij" = Eij, "processedData" = small_tidy, "drugList" = drug_keep, "AEList" = AE_keep, "N_E_list" = N_E_list)

  } else if (data.squashing == TRUE) {

    #if (!require('openEBGM')) install.packages('openEBGM'); library('openEBGM')

    small_rename = small_tidy[, c("primaryid", "PT", "prod_ai")]
    colnames(small_rename) = c("id", "var1", "var2")

    small_rename = subset(small_rename, var1 %in% AE_keep & var2 %in% drug_keep)
    small_rename = na.omit(small_rename)
    small_rename = droplevels(small_rename)

    Nij_s = Nij[which(Nij$frequency != 0), ]
    Eij_s = Eij[which(Nij$frequency != 0), ]

    Eij_prime = Eij_s[order(Eij_s$PT),]

    small_E_RR_PRR = openEBGM::processRaw(small_rename, zeroes = zeroes)
    if (zeroes == FALSE){
      small_E_RR_PRR$E = Eij_prime$baseline
    } else {
      small_E_RR_PRR$E = Eij$baseline
    }
    squashed = openEBGM::autoSquash(small_E_RR_PRR, keep_pts = keep_pts)

    drug_subsets <- split(small_E_RR_PRR, small_E_RR_PRR$var2)
    N_ij_squashed = list()

    for (i in names(drug_subsets)){
      N_ij_squashed[[i]] = openEBGM::autoSquash(drug_subsets[[i]], keep_pts = 30)
    }

    N_E_list = list()
    for (j in (1:ncol(N_ij))){
      jj = as.data.frame(cbind(N_ij[, j], E_ij[, j], rep(1, nrow(N_ij))))
      colnames(jj) = c("N", "E", "weight")
      N_E_list[[j]] = jj
    }

    NnE <- list("N_ij" = N_ij, "Nij" = Nij, "E_ij" = E_ij, "Eij" = Eij, "processedData" = small_tidy, "drugList" = drug_keep, "AEList" = AE_keep, "N_ij_squashed" = N_ij_squashed, "N_E_list" = N_E_list, "squashed" = squashed)


  }
  return(NnE)

}





