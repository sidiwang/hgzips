#' HGZIPS - MGPS posterior lambda
#'
#' This HZINB function.........
#'
#' @import stats
#' @param alpha1 the final estimation of alpha1 in ZINB
#' @param beta1 the final estimation of beta1 in ZINB
#' @param alpha2 the final estimation of alpha2 in ZINB
#' @param beta2 the final estimation of beta2 in ZINB
#' @param pi the final estimation of pi in ZINB
#' @param N vector of Nij values
#' @param E vector of Eij values
#'
#' @return the posterior mean of the lambda in HZINB
#' @export
#' @seealso
#'
#'
###########################################################
## Posterior expectation of lambda|N (MGPS)
###########################################################

# input N and E format: Nij$frequency, Eij$baseline

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

Posteror_MGPS_log = function(alpha1, beta1, alpha2, beta2, pi, N, E){
  # calculate Qx: the posterior probability that lambda came fromt he first component of the mixtrue, given N = x.
  Qx_1 = pi*dnbinom(N, size = alpha1, prob = beta1/(E + beta1))
  Qx_2 = (1 - pi)*dnbinom(N, size = alpha2, prob = beta2/(E + beta2))
  Qx = Qx_1/(Qx_1 + Qx_2)

  # calculate posterior expectation of lambda|N:
  post_mean_lambda = Qx*((digamma(alpha1 + N) - log(beta1) + E)) + (1 - Qx)*((digamma(alpha2 + N) - log(beta2) + E))

  return(post_mean_lambda)

}
