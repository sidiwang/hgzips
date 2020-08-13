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

# input, N = Nij$frequency, E = Eij$baseline

post_mean_lambda_ZINB = function(alpha, beta, N, E){
  post_mean_lambda = (alpha + Nij$frequency)/(beta + Eij$baseline)
  return(post_mean_lambda)
}

post_mean_loglambda_ZINB = function(alpha, beta, N, E){
  post_mean_loglambda = digamma(alpha + Nij$frequency) - log(beta + Eij$baseline)
  return(post_mean_loglambda)
}
