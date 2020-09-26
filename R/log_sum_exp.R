#' log_sum_exp and LSE_R
#'
#' This \code{log_sum_exp} function calculates the log sum of exponentials of two numbers. This trick is to increase accuracy and avoid underflow and overflow problems when very small or very large numbers are represented directly.
#'
#' @name log_sum_exp
#' @aliases log_sum_exp
#' @param u first number
#' @param v second number
#' @rdname log_sum_exp
#' @aliases log_sum_exp
#' @return \code{log_sum_exp} returns logged sum of exponential of u and v
#' @examples
#' log_sum_exp(1, 2)
#' LSE_R(c(1,2,3,4))
#' @export
log_sum_exp = function(u, v) {
  sum =  max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))
  return(sum)
}



#' log Sum Exp of vector
#' @rdname log_sum_exp
#' @aliases LSE_R
#' @param vec vector to be summed
#' @return \code{LSE_R} the logged sum of the exponentials of all values in vector

#' @export

LSE_R <- function(vec){
  n.vec <- length(vec)
  vec <- sort(vec, decreasing = TRUE)
  Lk <- vec[1]
  for (k in 1:(n.vec-1)) {
    Lk <- max(vec[k+1], Lk) + log1p(exp(-abs(vec[k+1] - Lk)))
  }
  return(Lk)
}
