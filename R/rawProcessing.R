############################################################
## raw data processing, generating N_ij, E_ij, Nij, Eij
############################################################

# rawdata has four columns: primaryid, PT, drugname, prod_ai
# input: raw data
# output: two formats of Nij, two formats of Eij


rawProcessing = function(rawdata, frequency_threshold){
#  AE20Q1small <- rawdata
#  AE20Q1small$PT = tolower(AE20Q1small$PT)
#  AE20Q1small$prod_ai = tolower(AE20Q1small$prod_ai)
  
#  library(tidyr)
#  AE20Q1small_tidy = tidyr::separate_rows(AE20Q1small, prod_ai, sep = "\\\\")
  
  AE20Q1small_tidy = unique(rawdata[, c("primaryid", "PT", "prod_ai")])
  
  # list of unique PT
  PT_list = unique(AE20Q1small_tidy$PT)
  # list of unique drug
  prod_list = sort(unique(AE20Q1small_tidy$prod_ai))
  
  # raw N_ij matrix
  N_ij_raw = table(AE20Q1small_tidy$PT, AE20Q1small_tidy$prod_ai)
  
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
  
  NnE <- list("N_ij" = N_ij, "Nij" = Nij, "E_ij" = E_ij, "Eij" = Eij, "processedData" = AE20Q1small_tidy, "drugList" = drug_keep, "AEList" = AE_keep)
  return(NnE)
  
}
