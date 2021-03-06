#' HGZIPS - HZINB (assuming independence)
#'
#' This \code{HZINB_independence} function finds hyperparameter estimates by implementing the Expectation-Maximization (EM) algorithm and hierarchical zero-inflated negative binomial model with one gamma component.
#'
#' @name HZINB_independence
#' @import pscl
#' @import stats
#' @import emdbook
#'
#' @param grid_a alpha value grid
#' @param grid_b beta value grid
#' @param grid_omega omega value grid
#' @param init_pi_k initial probability of each alpha value for implementing the EM algorithm
#' @param init_pi_l inital probability of each beta value for implementing the EM algorithm
#' @param init_pi_h initial probability of each omega value for implementing the EM algprithm
#' @param dataset a list of squashed datasets that include N_ij, E_ij and weights for each drug (j). This dataset list can be generated by the rawProcessing function in this package.
#' @param iteration number of EM algorithm iterations to run
#' @param Loglik whether to return the loglikelihood of each iteration or not (TRUE or FALSE)
#' @param zeroes A logical scalar specifying if zero counts should be included.
#' @param N_star the minimum Nij count size to be used for hyperparameter estimation. If zeroes are included in Nij vector, please set N_star = NULL
#'


# +-x +-x +-x +-x +-x +-x +-x +-x
#  assuming independence
# +-x +-x +-x +-x +-x +-x +-x +-x

#' HZINB_independence
#' @return \code{HZINB_independence} a list of estimated probability of each alpha, beta, omega combination and their corresponding loglikelihood (optional)
#' \itemize{
#' \item \code{theta_EM} Estimate of hyperparameters for each EM iteration
#' \item \code{llh} logliklihood for each EM iteration (optional)
#' }
#' @export
HZINB_independence = function(grid_a, grid_b, grid_omega, init_pi_k, init_pi_l, init_pi_h, dataset, iteration, Loglik = FALSE, zeroes = FALSE, N_star = 1){
  ## EM algorithm

  for (k in 1:length(dataset)){
    if (!is.null(N_star)){
      dataset[[k]] = subset(dataset[[k]], N >= N_star)
    }
  }

  if (zeroes == FALSE){
    K = length(grid_a)
    L = length(grid_b)
    #grid_omega = grid_omega

    #  if (!require('countreg')) install.packages('countreg'); library('countreg')

    all_combinations = as.data.frame(matrix(NA, K*L, 2))
    colnames(all_combinations) = c("a_j", "b_j")
    all_combinations$a_j = rep(grid_a, 10)

    for (i in c(1:L)){
      all_combinations$b_j[((i - 1)*10 + 1):(i*10)] = rep(grid_b[i], 10)
    }

    ## EM algorithm

    # initialization
    N.EM <- iteration  # number of E-M iterations
    #iteration_50 = pi_klh[50,]

    pi_klh_K = matrix(NA, N.EM + 2, K)
    pi_klh_L = matrix(NA, N.EM + 2, L)

    pi_klh_K[1,] = init_pi_k
    pi_klh_L[1,] = init_pi_l

    pi_klh_all_combinations = as.data.frame(matrix(NA, K*L, 2))
    colnames(pi_klh_all_combinations) = c("K", "L")

    denominator = rep(NA, length(dataset))
    numerator = rep(NA, length(dataset))
    #ratio = rep(NA, ncol(N_ij))
    joint_probs = as.data.frame(matrix(NA, nrow(all_combinations), length(dataset)))

    for (j in 1:length(dataset)){
      for (m in 1:nrow(all_combinations)){
        joint_probs[m,j] = sum(dataset[[j]]$weight * (dnbinom(dataset[[j]]$N, size = all_combinations$a_j[m], prob = all_combinations$b_j[m]/(dataset[[j]]$E + all_combinations$b_j[m]), log = TRUE) - log1p(-pnbinom(N_star - 1, size = all_combinations$a_j[m], prob = all_combinations$b_j[m]/(dataset[[j]]$E + all_combinations$b_j[m]), log = TRUE))))
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

    ratio = as.data.frame(matrix(NA, K*L, length(dataset)))

    for (i in 1:(N.EM + 1)) {

      pi_klh_all_combinations$K = rep(pi_klh_K[i,], 10)

      for (ii in c(1:L)){
        pi_klh_all_combinations$L[((ii - 1)*10 + 1):(ii*10)] = rep(pi_klh_L[i,][ii], 10)
      }

      pi_klh_all_combinations$prod = pi_klh_all_combinations$K * pi_klh_all_combinations$L

      for (m in 1:nrow(all_combinations)){
        for (j in 1:length(dataset)){
          denominator[j] = LSE_R(log(pi_klh_all_combinations$prod) + joint_probs[,j])
          numerator[j] = log(pi_klh_K[i, which(grid_a == all_combinations$a_j[m])]*pi_klh_L[i, which(grid_b == all_combinations$b_j[m])]) + joint_probs[m, j]
          ratio[m,j] = numerator[j] - denominator[j]
        }
      }

      # if (Loglik == TRUE){
      #    RATIO[[i]] = ratio
      #  } else {
      #    RATIO = NULL
      #  }

      all = cbind(all_combinations, ratio)

      for (iv in 1:nrow(ratio)){
        all$Sum[iv] = LSE_R(ratio[iv,])
      }

      all$Sum = unlist(all$Sum)
      temp = subset(all, !is.na(Sum))
      overallSum = LSE_R(temp$Sum)
      sum_a_j = aggregate(temp$Sum, by = list(Category = temp$a_j), FUN=LSE_R)
      sum_b_j = aggregate(temp$Sum, by = list(Category = temp$b_j), FUN=LSE_R)

      a_id = NULL
      b_id = NULL
      omega_id = NULL

      for (kk in 1:length(grid_a)){
        a_id = append(a_id, ifelse(sum(grid_a[kk] == sum_a_j$Category) == 0, kk, next))
      }

      for (kk in 1:length(grid_b)){
        b_id = append(b_id, ifelse(sum(grid_b[kk] == sum_b_j$Category) == 0, kk, next))
      }

      if (length(a_id) == 0){
        pi_klh_K[i + 1, ] = exp(sum_a_j$x - overallSum)
      } else {
        pi_klh_K[i + 1, ][-a_id] = exp(sum_a_j$x - overallSum)
        pi_klh_K[i + 1, ][a_id] = 0
      }

      if (length(b_id) == 0){
        pi_klh_L[i + 1, ] = exp(sum_b_j$x - overallSum)
      } else {
        pi_klh_L[i + 1, ][-b_id] = exp(sum_b_j$x - overallSum)
        pi_klh_L[i + 1, ][b_id] = 0
      }

    }

    result = list("pi_K" = pi_klh_K[-(N.EM + 2), ], "pi_L" = pi_klh_L[-(N.EM + 2), ])

  } else {

  K = length(grid_a)
  L = length(grid_b)
  H = length(grid_omega)

  #install.packages("countreg", repos="http://R-Forge.R-project.org")
  #library(countreg)

  all_combinations = as.data.frame(matrix(NA, K*L*H, 3))
  colnames(all_combinations) = c("a_j", "b_j", "omega_j")
  all_combinations$a_j = rep(grid_a, 100)

  for (i in c(1:L)){
    all_combinations$b_j[((i - 1)*100 + 1):(i*100)] = rep(grid_b[i], 100)
  }

  for (i in c(1:H)){
    row_num = c(((i - 1)*10 + 1):(i*10))
    all_combinations$omega_j[row_num] = rep(grid_omega[i], 10)
  }

  all_combinations$omega_j = rep(all_combinations$omega_j[1:100], 10)

  # initialization
  N.EM <- iteration  # number of E-M iterations

  pi_klh_K = matrix(NA, N.EM + 2, K)
  pi_klh_L = matrix(NA, N.EM + 2, L)
  pi_klh_H = matrix(NA, N.EM + 2, H)

  pi_klh_K[1,] = init_pi_k
  pi_klh_L[1,] = init_pi_l
  pi_klh_H[1,] = init_pi_h


  pi_klh_all_combinations = as.data.frame(matrix(NA, K*L*H, 3))
  colnames(pi_klh_all_combinations) = c("K", "L", "H")

  denominator = rep(NA, length(dataset))
  numerator = rep(NA, length(dataset))
  ratio = as.data.frame(matrix(NA, K*H*L, length(dataset)))
  #llh_col = rep(NA, ncol(N_ij))
  #RATIO = list()

  joint_probs = as.data.frame(matrix(NA, nrow(all_combinations), length(dataset)))

  for (j in 1:length(dataset)){
    for (m in 1:nrow(all_combinations)){
      joint_probs[m,j] = sum(dataset[[j]]$weight * emdbook::dzinbinom(dataset[[j]]$N, mu = (dataset[[j]]$E/all_combinations$b_j[m])*all_combinations$a_j[m], size = all_combinations$a_j[m], zprob = all_combinations$omega_j[m], log = TRUE))
    }
  }

  llh_j = rep(NA, length(dataset))
  llh = rep(NA, N.EM + 1)

  for (i in 1:(N.EM + 1)) {

    pi_klh_all_combinations$K = rep(pi_klh_K[i,], 100)

    for (ii in c(1:L)){
      pi_klh_all_combinations$L[((ii - 1)*100 + 1):(ii*100)] = rep(pi_klh_L[i,][ii], 100)
    }

    for (iii in c(1:H)){
      row_num = c(((iii - 1)*10 + 1):(iii*10))
      pi_klh_all_combinations$H[row_num] = rep(pi_klh_H[i,][iii], 10)
    }

    pi_klh_all_combinations$H = rep(pi_klh_all_combinations$H[1:100], 10)
    pi_klh_all_combinations$prod = pi_klh_all_combinations$K * pi_klh_all_combinations$L * pi_klh_all_combinations$H

    for (m in 1:nrow(all_combinations)){
      for (j in 1:length(dataset)){
        denominator[j] = LSE_R(log(pi_klh_all_combinations$prod) + joint_probs[,j])
        numerator[j] = log(pi_klh_K[i, which(grid_a == all_combinations$a_j[m])]*pi_klh_L[i, which(grid_b == all_combinations$b_j[m])]*pi_klh_H[i, which(grid_omega == all_combinations$omega_j[m])]) + joint_probs[m, j]
        ratio[m,j] = numerator[j] - denominator[j]
      }
    }

   # if (Loglik == TRUE){
  #    RATIO[[i]] = ratio
  #  } else {
  #    RATIO = NULL
  #  }

    all = cbind(all_combinations, ratio)

    for (iv in 1:nrow(ratio)){
      all$Sum[iv] = LSE_R(ratio[iv,])
    }

    all$Sum = unlist(all$Sum)
    temp = subset(all, !is.na(Sum))
    overallSum = LSE_R(temp$Sum)
    sum_a_j = aggregate(temp$Sum, by = list(Category = temp$a_j), FUN=LSE_R)
    sum_b_j = aggregate(temp$Sum, by = list(Category = temp$b_j), FUN=LSE_R)
    sum_omega_j = aggregate(temp$Sum, by = list(Category = temp$omega_j), FUN=LSE_R)

    a_id = NULL
    b_id = NULL
    omega_id = NULL

    for (kk in 1:length(grid_a)){
      a_id = append(a_id, ifelse(sum(grid_a[kk] == sum_a_j$Category) == 0, kk, next))
    }

    for (kk in 1:length(grid_b)){
      b_id = append(b_id, ifelse(sum(grid_b[kk] == sum_b_j$Category) == 0, kk, next))
    }

    for (kk in 1:length(grid_omega)){
      omega_id = append(omega_id, ifelse(sum(grid_omega[kk] == sum_omega_j$Category) == 0, kk, next))
    }

    if (length(a_id) == 0){
      pi_klh_K[i + 1, ] = exp(sum_a_j$x - overallSum)
    } else {
      pi_klh_K[i + 1, ][-a_id] = exp(sum_a_j$x - overallSum)
      pi_klh_K[i + 1, ][a_id] = 0
    }

    if (length(b_id) == 0){
      pi_klh_L[i + 1, ] = exp(sum_b_j$x - overallSum)
    } else {
      pi_klh_L[i + 1, ][-b_id] = exp(sum_b_j$x - overallSum)
      pi_klh_L[i + 1, ][b_id] = 0
    }

    if (length(omega_id) == 0){
      pi_klh_H[i + 1, ] = exp(sum_omega_j$x - overallSum)
    } else {
      pi_klh_H[i + 1, ][-omega_id] = exp(sum_omega_j$x - overallSum)
      pi_klh_H[i + 1, ][omega_id] = 0
    }



    if (Loglik == TRUE){

      for (j in 1:length(dataset)){

        pi_klh_all_combinations$logSum = log(pi_klh_all_combinations$K) + log(unlist(pi_klh_all_combinations$L)) + log(unlist(pi_klh_all_combinations$H)) + joint_probs[,j]
        llh_j[j] = LSE_R(pi_klh_all_combinations$logSum)

      }

      llh[i] = sum(llh_j)
      print(i)

    } else {
      llh = NULL
    }
  }

  result = list("pi_K" = pi_klh_K[-(N.EM + 2), ], "pi_L" = pi_klh_L[-(N.EM + 2), ], "pi_H" = pi_klh_H[-(N.EM + 2), ], "Loglik" = llh)
  }
  return(result)

}






