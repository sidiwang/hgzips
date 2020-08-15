#'
#' log Sum Exp of vector
#'
#' @param vec vector to be summed
#' @return the logged sum of the exponentials of all values in vector
#' @export
#' @seealso
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
