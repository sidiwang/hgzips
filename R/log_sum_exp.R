#'
#'log Sum Exp
#' @name LogSumExp
#'




#' @param u first number
#' @param v second number
#' @rdname LogSumExp
#' @aliases log_sum_exp
#' @return logged sum of exponential of u and v
#' @export
log_sum_exp = function(u, v) {
  sum =  max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))
  return(sum)
}



#' log Sum Exp of vector
#' @rdname LogSumExp
#' @aliases LSE_R
#' @param vec vector to be summed
#' @return the logged sum of the exponentials of all values in vector
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
