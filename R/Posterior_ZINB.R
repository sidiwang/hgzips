#' HGZIPS - ZINB posterior lambda
#'
#' This HZINB function.........
#' @name ZINB_posterior
#' @import stats
#'
#' @param alpha the final estimation of alpha in ZINB
#' @param beta the final estimation of beta in ZINB
#' @param N vector of Nij values
#' @param E vector of Eij values
#' @return the posterior mean of the lambda in ZINB
#' @seealso
#'
###########################################################
## Posterior expectation of lambda|N (MGPS)
###########################################################

# input, N = Nij$frequency, E = Eij$baseline
#' @rdname ZINB_posterior
#' @aliases post_mean_lambda_ZINB
#' post_mean_lambda_ZINB
#' @export
post_mean_lambda_ZINB = function(alpha, beta, N, E){
  post_mean_lambda = (alpha + N)/(beta + E)
  return(post_mean_lambda)
}

#' @rdname ZINB_posterior
#' @aliases post_mean_loglambda_ZINB
#' post_mean_loglambda_ZINB
#' @return posterior mean of logged lambda
#' @export
post_mean_loglambda_ZINB = function(alpha, beta, N, E){
  post_mean_loglambda = digamma(alpha + N) - log(beta + E)
  return(post_mean_loglambda)
}
