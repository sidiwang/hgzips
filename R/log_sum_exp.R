
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




log_sum_exp = function(u, v) {
  sum =  max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))
  return(sum)
}
