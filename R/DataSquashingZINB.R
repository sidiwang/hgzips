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


################################################
##  data squashing ZINB one gamma
################################################


LSE_R <- function(vec){
  n.vec <- length(vec)
  vec <- sort(vec, decreasing = TRUE)
  Lk <- vec[1]
  for (k in 1:(n.vec-1)) {
    Lk <- max(vec[k+1], Lk) + log1p(exp(-abs(vec[k+1] - Lk)))
  }
  return(Lk)
}

DataSquashingZINB = function(alpha, beta, omega, squashedData, iteration, Loglik){

  N_E_weight = squashedData
  N_E_weight$logisticY = ifelse(N_E_weight$N == 0, 0, 1)


  # initialization
  N.EM <- iteration  # number of E-M iterations
  theta_EM = matrix(rep(0, 3 * (N.EM + 1)), nrow = N.EM + 1, ncol = 3)  # matrix for holding theta in the E-M iterations
  theta_EM = as.data.frame(theta_EM)
  theta_EM[1, ] <- c(alpha, beta, omega)  # initial value for parameters
  colnames(theta_EM) = c("alpha", "beta", "omega")

  T_ij = rep(NA, nrow(N_E_weight))
  cluster = rep(NA, nrow(N_E_weight))
  #u_ij = rep(NA, nrow(Nij))


  # count = 0
  zero = which(N_E_weight$N == 0)
  #u_ij = ifelse(Nij$frequency == 0, 1, 0)

  # the E-M iterations
  for (i in c(1:N.EM)) {
    cat("Iteration", i, "\n")

    # Expectation step:
    T_ij_zero_1 = log(theta_EM$omega[i])
    T_ij_zero_2 = log(theta_EM$omega[i] + (1 - theta_EM$omega[i])*(((theta_EM$beta[i])/(theta_EM$beta[i] + N_E_weight$E[zero]))^(theta_EM$alpha[i])))
    T_ij[zero] = exp(T_ij_zero_1 - T_ij_zero_2)
    T_ij[-zero] = 0

    # Maximization step:
    theta_EM$omega[i + 1] = sum(T_ij*N_E_weight$weight)/(sum(T_ij*N_E_weight$weight) + sum((1 - T_ij)*N_E_weight$weight))

    ZNegBin = function(parameter, N, E, Ts){
      alpha = parameter[1]
      beta = parameter[2]

      max = sum((1 - Ts)*N_E_weight$weight*dnbinom(N, size = alpha, prob = beta/(E + beta), log = TRUE))

      return(-max)
    }

    theta_EM[i + 1, 1:2] = optim(c(theta_EM$alpha[i], theta_EM$beta[i]), ZNegBin, method = "BFGS", N = N_E_weight$N, E = N_E_weight$E, Ts = T_ij)$par

  }


  if (Loglik == TRUE){

    llh1 = rep(NA, nrow(theta_EM))
    log.dens = rep(0, length(N_E_weight$N))

    for (i in c(1: nrow(theta_EM))){

      log.dens = ifelse(N == 0, log(theta_EM$omega[i] + exp(log(1 - theta_EM$omega[i]) + dnbinom(0, size = theta_EM$alpha[i], prob = theta_EM$beta[i]/(E + theta_EM$beta[i]), log = TRUE))), log(1 - theta_EM$omega[i]) + dnbinom(N, size = theta_EM$alpha[i], prob = theta_EM$beta[i]/(E + theta_EM$beta[i]), log = TRUE))

      llh1[i] = sum(log.dens*N_E_weight$weight)
      }
    } else {
    llh1 = NULL
  }

  result = list("theta_EM" = theta_EM, "Loglik" = llh1)
  return(result)

}




