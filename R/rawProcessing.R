#' HGZIPS - rawProcessing
#'
#' This function.........
#'
#' @import stats
#' @import openEBGM
#'
#' @param rawdata rawdata input
#' @param frequency_threshold minimum......
#' @param data.squashing whether to conduct data squashing or not (TRUE or FALSE)

#' @return a list of objects after raw data processing
#' @export
#' @seealso

############################################################
## raw data processing, generating N_ij, E_ij, Nij, Eij
############################################################

# rawdata has four columns: primaryid, PT, drugname, prod_ai
# input: raw data
# output: two formats of Nij, two formats of Eij


rawProcessing = function(rawdata, frequency_threshold, data.squashing = FALSE){


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
  N_ij = as.data.frame.matrix(N_ij_raw[-rows_to_drop, -columns_to_drop])

  drug_keep = colnames(N_ij)
  AE_keep = rownames(N_ij)

  # another format
  Nij = as.data.frame(N_ij_raw[-rows_to_drop, -columns_to_drop])
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

  # another format
  Eij = as.data.frame(as.table(E_ij))
  colnames(Eij) = c("PT", "prod_ai", "baseline")

  E_ij = as.data.frame(E_ij)

  if (data.squashing == FALSE){

     NnE <- list("N_ij" = N_ij, "Nij" = Nij, "E_ij" = E_ij, "Eij" = Eij, "processedData" = small_tidy, "drugList" = drug_keep, "AEList" = AE_keep)

     } else if (data.squashing == TRUE) {

    if (!require('openEBGM')) install.packages('openEBGM'); library('openEBGM')

    small_rename = small_tidy[, c("primaryid", "PT", "prod_ai")]
    colnames(small_rename) = c("id", "var1", "var2")

    small_rename = subset(small_rename, var1 %in% AE_keep & var2 %in% drug_keep)
    small_rename = na.omit(small_rename)
    small_rename = droplevels(small_rename)

    Eij_prime = Eij[order(Eij$PT),]

    small_E_RR_PRR = openEBGM::processRaw(small_rename, zeroes = TRUE)
    small_E_RR_PRR$E = Eij_prime$baseline
    squashed = openEBGM::autoSquash(small_E_RR_PRR, keep_pts = 190000)

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





