#' HGZIPS - HZINB posterior lambda
#'
#' This HZINB function.........
#'
#' @import stats
#'
#' @param posterior_a_j the posterior average of HZINB a_j, can be produced by the posterior_abomega_HZINB function in the package
#' @param posterior_b_j the posterior average of HZINB b_j, can be produced by the posterior_abomega_HZINB function in the package
#' @param posterior_omega_j the posterior average of HZINB b_j, can be produced by the posterior_abomega_HZINB function in the package
#' @param N_ij matrix of N_ij, i = AE, j = drugs
#' @param E_ij matrix of E_ij, i = AE, j = drugs

#' @return the posterior mean of the lambda in HZINB
#' @export
#' @seealso
#'
##########################################################
## Posterior HZINB One Gamma
##########################################################

# N = N_ij, E = E_ij

posteriorHZINB = function(posterior_a_j, posterior_b_j, posterior_omega_j, N_ij, E_ij){

  post_HZINB = matrix(NA, nrow(N_ij), ncol(N_ij))

  for (j in 1 : ncol(N_ij)){
    post_HZINB[, j] =  (posterior_a_j[j] + N_ij[, j])/(posterior_b_j[j] + E_ij[, j])
    Qx = posterior_omega_j[j]/(posterior_omega_j[j] + (1 - posterior_omega_j[j])*dnbinom(0, size = posterior_a_j[j], prob = posterior_b_j[j]/E_ij[, j] + posterior_b_j[j]))
    post_HZINB[which(N_ij[, j] == 0, arr.ind = TRUE)] = (Qx*posterior_a_j[j]/posterior_b_j[j] + (1 - Qx)*posterior_a_j[j]/(posterior_b_j[j] + E_ij[, j]))[which(N_ij[, j] == 0, arr.ind = TRUE)]

  }

  return(post_HZINB)
}



