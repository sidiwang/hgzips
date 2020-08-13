#' Create a complete ggplot appropriate to a particular data type
#'
#' \code{autoplot} uses ggplot2 to draw a particular plot for an object of a particular class in a single command.
#' This defines the S3 generic that other classes and packages can extend.
#'
#' @param object an object, whose class will determin the behaviour of autoplot
#' @param  ...other auguments passed to specific methods
#' @return a ggplot object
#' @export
#' @seealso  \code{\link{ggplot}} and \code{\link{fortify}}
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
