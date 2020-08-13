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






LSE_R <- function(vec){
  n.vec <- length(vec)
  vec <- sort(vec, decreasing = TRUE)
  Lk <- vec[1]
  for (k in 1:(n.vec-1)) {
    Lk <- max(vec[k+1], Lk) + log1p(exp(-abs(vec[k+1] - Lk)))
  }
  return(Lk)
}
