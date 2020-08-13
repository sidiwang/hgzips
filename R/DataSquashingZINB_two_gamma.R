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



# input N = Nij$frequency, E = Eij$baseline


DataSquashingZINB_two_gamma = function(alpha1, beta1, alpha2, beta2, pi, omega, N, E, weight, iteration, Loglik){

  log_sum_exp = function(u, v) {
    sum =  max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))
    return(sum)
  }

  logisticY = rep(NA, length(N))
  # estimate omega_ij (logistic regression)
  logisticY = ifelse(N == 0, 0, 1)

  # initialization
  N.EM <- iteration  # number of E-M iterations
  theta_EM = matrix(rep(0, 6 * (N.EM + 1)), nrow = N.EM + 1, ncol = 6)  # matrix for holding theta in the E-M iterations
  theta_EM = as.data.frame(theta_EM)
  theta_EM[1, ] <- c(alpha1, beta1, alpha2, beta2, pi, omega)  # initial value for parameters
  colnames(theta_EM) = c("alpha1", "beta1", "alpha2", "beta2", "pi", "omega")

  T_ij = rep(NA, nrow(Nij))
  cluster = rep(NA, nrow(Nij))
  u_ij = rep(NA, nrow(Nij))


  # count = 0
  zero = which(N == 0)
  u_ij = ifelse(N == 0, 1, 0)

  # the E-M iterations
  for (i in c(1:N.EM)) {

    print(i)

    # Expectation step:
    T_ij_zero_1 = theta_EM$pi[i]*exp(dnbinom(0, size = theta_EM$alpha1[i], prob = theta_EM$beta1[i]/(E[zero] + theta_EM$beta1[i]), log = TRUE))
    T_ij_zero_2 = theta_EM$pi[i]*exp(dnbinom(0, size = theta_EM$alpha2[i], prob = theta_EM$beta2[i]/(E[zero] + theta_EM$beta2[i]), log = TRUE))
    T_ij[zero] = exp(log(theta_EM$omega[i] + (1 - theta_EM$omega[i])*T_ij_zero_1) - log(2*theta_EM$omega[i] + (1 - theta_EM$omega[i])*(T_ij_zero_1 + T_ij_zero_2)))

    probs1 <- theta_EM$beta1[i]/(E[-zero] + theta_EM$beta1[i])
    probs2 <- theta_EM$beta2[i]/(E[-zero] + theta_EM$beta2[i])
    T_ij_1 = dnbinom(N[-zero], size = theta_EM$alpha1[i], prob = probs1, log = TRUE)
    T_ij_2 = dnbinom(N[-zero], size = theta_EM$alpha2[i], prob = probs2, log = TRUE)

    logBF <- T_ij_2 - T_ij_1 + log(1 - theta_EM$pi[i]) - log(theta_EM$pi[i])
    T_ij[-zero][logBF < 0] <- exp(-log1p(exp(logBF[logBF < 0])))
    T_ij[-zero][logBF >= 0] <- exp(-logBF[logBF >= 0] - log1p(exp(-logBF[logBF>=0])))

    ZNegBin = function(parameter, u, N, E, weight, Ts){
      alpha1 = parameter[1]
      beta1 = parameter[2]
      alpha2 = parameter[3]
      beta2 = parameter[4]
      pi = parameter[5]
      omega = parameter[6]

      max = sum(Ts *weight*(u*(log(omega + (1 - omega)*(pi*exp(dnbinom(0, size = alpha1, prob = beta1/(E + beta1), log = TRUE))))) + (1 - u)*(log(pi) + log(1 - omega) + dnbinom(N, size = alpha1, prob = beta1/(E + beta1), log = TRUE)))) +
        sum((1 - Ts) *weight*(u*(log(omega + (1 - omega)*((1 - pi)*exp(dnbinom(0, size = alpha2, prob = beta2/(E + beta2), log = TRUE))))) + (1 - u)*(log(1 - pi) + log(1 - omega) + dnbinom(N, size = alpha2, prob = beta2/(E + beta2), log = TRUE))))

      return(-max)
    }

    theta_EM[i + 1, 1:6] = optim(c(theta_EM$alpha1[i], theta_EM$beta1[i], theta_EM$alpha2[i], theta_EM$beta2[i], theta_EM$pi[i], theta_EM$omega[i]), ZNegBin, method = "BFGS", u = u_ij, N = N, E = E, weight = weight, Ts = T_ij)$par

  }

  if (Loglik == TRUE){

    llh1 = rep(NA, nrow(theta_EM))
    log.dens = rep(0, length(N))

    for (i in 1:nrow(theta_EM)){
      #probs1 <- theta_EM$beta1[i]/(E[-zero] + theta_EM$beta1[i])
      #probs2 <- theta_EM$beta2[i]/(E[-zero] + theta_EM$beta2[i])

      log.dens = ifelse(N == 0, log_sum_exp(log(theta_EM$omega[i]), (log(1 - theta_EM$omega[i]) + log_sum_exp((log(theta_EM$pi[i]) + dnbinom(0, size = theta_EM$alpha1[i], prob = theta_EM$beta1[i]/(E + theta_EM$beta1[i]), log = TRUE)), (log(1 - theta_EM$pi[i]) + dnbinom(0, size = theta_EM$alpha2[i], prob = theta_EM$beta2[i]/(E + theta_EM$beta2[i]), log = TRUE))))), (log(1 - theta_EM$omega[i]) + log_sum_exp((log(theta_EM$pi[i]) + dnbinom(N, size = theta_EM$alpha1[i], prob = theta_EM$beta1[i]/(E + theta_EM$beta1[i]), log = TRUE)), (log(1 - theta_EM$pi[i]) + dnbinom(N, size = theta_EM$alpha2[i], prob = theta_EM$beta2[i]/(E + theta_EM$beta2[i]), log = TRUE)))))

      llh1[i] = sum(log.dens*weight)

      }
  } else {
    llh1 = NULL
  }

  result = list("theta_EM" = theta_EM, "Loglik" = llh1)
  return(result)

}





