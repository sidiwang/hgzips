################################################
##  data squashing ZINB two gamma
################################################

if (!require('openEBGM')) install.packages('openEBGM'); library('openEBGM')

AE20Q1small_rename = AE20Q1small_tidy[, c("primaryid", "PT", "prod_ai")]
colnames(AE20Q1small_rename) = c("id", "var1", "var2")

AE20Q1small_rename = subset(AE20Q1small_rename, var1 %in% AE_keep & var2 %in% drug_keep)
AE20Q1small_rename = na.omit(AE20Q1small_rename)
AE20Q1small_rename = droplevels(AE20Q1small_rename)

Eij_prime = Eij[order(Eij$PT),]

AE20Q1small_E_RR_PRR = processRaw(AE20Q1small_rename, zeroes = TRUE)
AE20Q1small_E_RR_PRR$E = Eij_prime$baseline
squashed = autoSquash(AE20Q1small_E_RR_PRR, keep_pts = 190000)

drug_subsets <- split(AE20Q1small_E_RR_PRR, AE20Q1small_E_RR_PRR$var2)
N_ij_squashed = list()

for (i in names(drug_subsets)){
  N_ij_squashed[[i]] = autoSquash(drug_subsets[[i]], keep_pts = 30)
}


#squashed = as.data.frame(cbind(AE20Q1small_E_RR_PRR$N, AE20Q1small_E_RR_PRR$E))
#colnames(squashed) = c("N", "E")
#squashed$weight = 1


##########################################################
## HZINB: check the range of parameters
##########################################################

# input N = N_ij, E = E_ij

parRangeCheck = function(N_ij, E_ij){

  if (!require('countreg')) install.packages('countreg'); library('countreg')

  a_j = rep(NA, ncol(N_ij))
  b_j = rep(NA, ncol(N_ij))
  omega_j = rep(NA, ncol(N_ij))

  for (j in 1:ncol(N_ij)){
    data = as.data.frame(cbind(N_ij[,j], rep(1, nrow(N_ij))))
    colnames(data) = c("N_ij", "X")
    tryCatch({
      model = zeroinfl(N_ij ~ 1, data = data, dist = "negbin", offset = log(E_ij[,j]))
      a_j[j] = model$theta
      omega_j[j] = exp(coef(model)[2])/(1 + exp(coef(model)[2]))
      b_j[j] = a_j[j]/exp(model$coefficients[1]$count)

    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }

  result = list("a_j" = a_j, "b_j" = b_j, "omega_j" = omega_j)
  return(result)

}



##########################################################
## HZINB: build a parameter grid
##########################################################

grid_HZINB = function(a_j, b_j, omega_j, K, L, H){

  grid = as.data.frame(matrix(NA, max(c(K, L, H)), 3))
  colnames(grid) = c("a_j", "b_j", "omega_j")

  for (i in c(1:H)){
    grid[i,3] = exp(log(omega_j[which.min(omega_j)]) + (i)/(K + 1)*(log(omega_j[which.max(omega_j)])  - log(omega_j[which.min(omega_j)])))
  }

  for (i in c(1:K)){
    grid[i,1] = exp(log(a_j[which.min(a_j)]) + (i - 1)/(K + 1)*(log(quantile(a_j, 0.99, names = FALSE, na.rm = TRUE)) - log(a_j[which.min(a_j)])))
  }

  for (i in c(1:L)){
    grid[i,2] = exp(log(b_j[which.min(b_j)]) + (i - 1)/(K + 1)*(log(quantile(b_j, 0.99, names = FALSE, na.rm = TRUE)) - log(b_j[which.min(b_j)])))
  }

  return(grid)

}


#####################################
## Hierarchical ZINB (one gamma)
#####################################
# +-x +-x +-x +-x +-x +-x +-x +-x
#  Not assuming independence
# +-x +-x +-x +-x +-x +-x +-x +-x


DataSquashingHZINB_one_gamma = function(K, L, H, grid, init_pi_klh, N_ij_squashed, iteration, Loglik){

  #install.packages("countreg", repos="http://R-Forge.R-project.org")
  library(countreg)

  all_combinations = as.data.frame(matrix(NA, K*L*H, 3))
  colnames(all_combinations) = c("a_j", "b_j", "omega_j")
  all_combinations$a_j = rep(grid$a_j, 100)

  for (i in c(1:L)){
    all_combinations$b_j[((i - 1)*100 + 1):(i*100)] = rep(grid$b_j[i], 100)
  }

  for (i in c(1:H)){
    row_num = c(((i - 1)*10 + 1):(i*10))
    all_combinations$omega_j[row_num] = rep(grid$omega_j[i], 10)
  }

  all_combinations$omega_j = rep(all_combinations$omega_j[1:100], 10)

  ## EM algorithm

  # initialization
  N.EM <- iteration  # number of E-M iterations
  #iteration_50 = pi_klh[50,]

  pi_klh = matrix(NA, N.EM + 1, K*L*H)
  pi_klh[1,] = init_pi_klh

  denominator = rep(NA, length(N_ij_squashed))
  numerator = rep(NA, length(N_ij_squashed))
  #ratio = rep(NA, ncol(N_ij))
  joint_probs = as.data.frame(matrix(NA, nrow(all_combinations), length(N_ij_squashed)))

  for (j in 1:length(N_ij_squashed)){
    for (m in 1:nrow(all_combinations)){
      joint_probs[m,j] = sum(N_ij_squashed[[j]]$weight * dzinbinom(N_ij_squashed[[j]]$N, mu = (N_ij_squashed[[j]]$E/all_combinations$b_j[m])*all_combinations$a_j[m], size = all_combinations$a_j[m], pi = all_combinations$omega_j[m], log = TRUE))
    }
  }

  LSE_R <- function(vec){
    n.vec <- length(vec)
    vec <- sort(vec, decreasing = TRUE)
    Lk <- vec[1]
    for (k in 1:(n.vec-1)) {
      Lk <- max(vec[k+1], Lk) + log1p(exp(-abs(vec[k+1] - Lk)))
    }
    return(Lk)
  }

  ratio = as.data.frame(matrix(NA, K*H*L, length(N_ij_squashed)))

  for (i in c(1 : N.EM)) {
    print(i)
    for (m in 1:nrow(all_combinations)){
      for (j in 1:length(N_ij_squashed)){
        ratio[m,j] = log(pi_klh[i, m]) + joint_probs[m,j] - LSE_R((log(pi_klh[i,]) + joint_probs[,j])[which((log(pi_klh[i,]) + joint_probs[,j]) != "NaN")])
      }
    }

    all = cbind(all_combinations, ratio)
    all$sum = rowSums(exp(ratio))
    overallSum = sum(all$sum, na.rm = TRUE)

    pi_klh[i + 1, ] = (all$sum)/overallSum
  }

  if (Loglik == TRUE){

    llh_j = rep(NA, ncol(N_ij))
    llh = rep(NA, N.EM)

    for (i in 1:N.EM){
      for (j in 1:ncol(N_ij)){
        llh_j[j] = unlist(LSE_R(log(pi_klh[i,]) + joint_probs[,j]))
      }
      llh[i] = sum(llh_j)
    }

  } else {
    llh = NULL
  }

  result = list("pi_klh" = pi_klh, "Loglik" = llh)
  return(result)

}



