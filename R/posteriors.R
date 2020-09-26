#' HGZIPS - calculating the posterior value of HZINB a_j, b_j and omega_j
#'
#' @name posterior
#'
#' @import stats
#' @import emdbook
#'
#' @param grid_a alpha value grid
#' @param grid_b beta value grid
#' @param grid_omega omega value grid
#' @param pi_klh_final_a_j final estimation of the probability of each alpha value, produced by the EM algorithm
#' @param pi_klh_final_b_j final estimation of the probability of each beta value, produced by the EM algorithm
#' @param pi_klh_final_omega_j final estimation of the probability of each omega value, produced by the EM algorithm
#' @param N_ij matrix of N_ij, i = AE, j = drugs
#' @param E_ij matrix of E_ij, i = AE, j = drugs

#' @seealso
#'
#'
#'
#'
#'



#' @rdname posterior
#' @return a list of estimated probability of each alpha, beta, omega combination and their corresponding loglikelihood (optional)
#' @export
#'
posterior_abol = function(grid_a, grid_b, grid_omega, pi_klh_final_a_j, pi_klh_final_b_j, pi_klh_final_omega_j, N_ij, E_ij){

  a_j = grid_a
  b_j = grid_b
  omega_j = grid_omega

  K = length(grid_a)
  L = length(grid_b)
  H = length(grid_omega)

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

  joint_probs = as.data.frame(matrix(NA, nrow(all_combinations), ncol(N_ij)))

  for (j in 1:ncol(N_ij)){
    for (m in 1:nrow(all_combinations)){
      joint_probs[m,j] = sum(emdbook::dzinbinom(N_ij[,j], mu = (E_ij[,j]/all_combinations$b_j[m])*all_combinations$a_j[m], size = all_combinations$a_j[m], zprob = all_combinations$omega_j[m], log = TRUE))
    }
  }

  pi_klh_all_combinations = as.data.frame(matrix(NA, K*L*H, 3))
  colnames(pi_klh_all_combinations) = c("K", "L", "H")
  pi_klh_all_combinations$K = rep(pi_klh_final_a_j, 100)

  for (i in c(1:L)){
    pi_klh_all_combinations$L[((i - 1)*100 + 1):(i*100)] = rep(pi_klh_final_b_j[i], 100)
  }

  for (i in c(1:H)){
    row_num = c(((i - 1)*10 + 1):(i*10))
    pi_klh_all_combinations$H[row_num] = rep(pi_klh_final_omega_j[i], 10)
  }

  pi_klh_all_combinations$H = rep(pi_klh_all_combinations$H[1:100], 10)

  pi_klh_all_combinations$prod = unlist(pi_klh_all_combinations$K) * unlist(pi_klh_all_combinations$L) * unlist(pi_klh_all_combinations$H)

  LSE_R <- function(vec){
    n.vec <- length(vec)
    vec <- sort(vec, decreasing = TRUE)
    Lk <- vec[1]
    for (k in 1:(n.vec-1)) {
      Lk <- max(vec[k+1], Lk) + log1p(exp(-abs(vec[k+1] - Lk)))
    }
    return(Lk)
  }

  expectation_a_j = rep(NA, ncol(N_ij))
  expectation_b_j = rep(NA, ncol(N_ij))
  expectation_omega_j = rep(NA, ncol(N_ij))
  expectation_lambda_j = matrix(NA, nrow(N_ij), ncol(N_ij))
  numerator_lambda_j = rep(NA, nrow(N_ij))

  for (j in 1 : ncol(N_ij)){

    print(j)

    denominator_a_j = LSE_R(joint_probs[, j] + log(pi_klh_all_combinations$prod))
    numerator_a_j = LSE_R(joint_probs[, j] + log(all_combinations$a_j) + log(pi_klh_all_combinations$prod))
    expectation_a_j[j] = exp(numerator_a_j - denominator_a_j)

    denominator_b_j = LSE_R(joint_probs[, j] + log(pi_klh_all_combinations$prod))
    numerator_b_j = LSE_R(joint_probs[, j] + log(all_combinations$b_j) + log(pi_klh_all_combinations$prod))
    expectation_b_j[j] = exp(numerator_b_j - denominator_b_j)

    denominator_omega_j = LSE_R(joint_probs[, j] + log(pi_klh_all_combinations$prod))
    numerator_omega_j = LSE_R(joint_probs[, j] + log(all_combinations$omega_j) + log(pi_klh_all_combinations$prod))
    expectation_omega_j[j] = exp(numerator_omega_j - denominator_omega_j)

    denominator_lambda_j = LSE_R(joint_probs[, j] + log(pi_klh_all_combinations$prod))
    for (i in 1 : nrow(N_ij)){
      numerator_lambda_j[i] = LSE_R(joint_probs[, j] + log((all_combinations$a_j + N_ij[i, j])/(all_combinations$b_j + E_ij[i, j])) + log(pi_klh_all_combinations$prod))
    }
    expectation_lambda_j[, j] = exp(numerator_lambda_j - denominator_lambda_j)

  }

  result = list("posterior_a_j" = expectation_a_j, "posterior_b_j" = expectation_b_j, "posterior_omega_j" = expectation_omega_j, "posterior_lambda_ij" = expectation_lambda_j)
  return(result)

}






# input, N = Nij$frequency, E = Eij$baseline
#' @rdname posterior
#' @aliases post_mean_lambda_ZINB
#' post_mean_lambda_ZINB
#' @export
post_mean_lambda_ZINB = function(alpha, beta, N, E){
  post_mean_lambda = (alpha + N)/(beta + E)
  return(post_mean_lambda)
}

#' @rdname posterior
#' @aliases post_mean_loglambda_ZINB
#' post_mean_loglambda_ZINB
#' @return posterior mean of logged lambda
#' @export
post_mean_loglambda_ZINB = function(alpha, beta, N, E){
  post_mean_loglambda = digamma(alpha + N) - log(beta + E)
  return(post_mean_loglambda)
}



#' @rdname posterior
#' @aliases Posteror_MGPS
#' @return a list of estimated probability of each alpha, beta, omega combination and their corresponding loglikelihood (optional)
#' @export
Posteror_MGPS = function(alpha1, beta1, alpha2, beta2, pi, N, E){
  # calculate Qx: the posterior probability that lambda came fromt he first component of the mixtrue, given N = x.
  Qx_1 = pi*dnbinom(N, size = alpha1, prob = beta1/(E + beta1))
  Qx_2 = (1 - pi)*dnbinom(N, size = alpha2, prob = beta2/(E + beta2))
  Qx = Qx_1/(Qx_1 + Qx_2)

  # calculate posterior expectation of lambda|N:
  post_mean_lambda = Qx*(alpha1 + N)/(beta1 + E) + (1 - Qx)*(alpha2 + N)/(beta2 + E)

  return(post_mean_lambda)

}

# input N and E format: Nij$frequency, Eij$baseline
#' @rdname posterior
#' @aliases Posteror_MGPS_log
#' @return a list of estimated probability of each alpha, beta, omega combination and their corresponding loglikelihood (optional)
#' @export
Posteror_MGPS_log = function(alpha1, beta1, alpha2, beta2, pi, N, E){
  # calculate Qx: the posterior probability that lambda came fromt he first component of the mixtrue, given N = x.
  Qx_1 = pi*dnbinom(N, size = alpha1, prob = beta1/(E + beta1))
  Qx_2 = (1 - pi)*dnbinom(N, size = alpha2, prob = beta2/(E + beta2))
  Qx = Qx_1/(Qx_1 + Qx_2)

  # calculate posterior expectation of lambda|N:
  post_mean_lambda = Qx*((digamma(alpha1 + N) - log(beta1) + E)) + (1 - Qx)*((digamma(alpha2 + N) - log(beta2) + E))

  return(post_mean_lambda)

}
