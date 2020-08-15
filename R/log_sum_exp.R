#'
#'log Sum Exp
#'
#'
#' @param u first number
#' @param v second number
#' @return logged sum of exponential of u and v
#' @export
#' @seealso
#'
#'




log_sum_exp = function(u, v) {
  sum =  max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))
  return(sum)
}
