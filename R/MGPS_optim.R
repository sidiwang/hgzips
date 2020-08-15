#' HGZIPS - MGPS (optim)
#'
#' This MGPS function.........
#'
#' @import stats
#'
#' @param alpha1 initial shape parameter value of the first gamma distribution for implementing the EM algprithm
#' @param beta1 initial rate parameter value of the first gamma distribution for implementing the EM algprithm
#' @param alpha2 initial shape parameter value of the second gamma distribution for implementing the EM algprithm
#' @param beta2 initial rate parameter value of the second gamma distribution for implementing the EM algprithm
#' @param pi initial xxxxxxx?
#' @param N vector of N_ij values
#' @param E vector of E_ij values
#' @param iteration number of EM algorithm iterations to run
#' @param Loglik whether to return the loglikelihood of each iteration or not (TRUE or FALSE)
#' @return a list of estimated parameters and their corresponding loglikelihood
#' @export
#' @seealso
#'
MGPS_optim = function(alpha1, beta1, alpha2, beta2, pi, N, E, iteration, Loglik){

  # initialization
  N.EM <- iteration  # number of E-M iterations
  theta_EM = matrix(rep(0, 5 * (N.EM + 1)), nrow = N.EM + 1, ncol = 5)  # matrix for holding theta in the E-M iterations
  theta_EM = as.data.frame(theta_EM)
  theta_EM[1, ] <- c(alpha1, beta1, alpha2, beta2, pi)  # initial value for parameters
  colnames(theta_EM) = c("alpha1", "beta1", "alpha2", "beta2", "pi")

  T_ij = rep(NA, nrow(Nij))
  cluster = rep(NA, nrow(Nij))

  # the E-M iterations
  for (i in c(1:N.EM)) {

    # Expectation step:
    probs1 <- theta_EM$beta1[i]/(E + theta_EM$beta1[i])
    probs2 <- theta_EM$beta2[i]/(E + theta_EM$beta2[i])
    T_ij_1 = dnbinom(N, size = theta_EM$alpha1[i], prob = probs1, log = TRUE)
    T_ij_2 = dnbinom(N, size = theta_EM$alpha2[i], prob = probs2, log = TRUE)

    logBF <- T_ij_2 - T_ij_1 + log(1 - theta_EM$pi[i]) - log(theta_EM$pi[i])
    T_ij[logBF < 0] <- exp(-log1p(exp(logBF[logBF < 0])))
    T_ij[logBF >= 0] <- exp(-logBF[logBF >= 0] - log1p(exp(-logBF[logBF>=0])))

    # Maximization step:

    theta_EM$pi[i + 1] = sum(T_ij)/length(N)

    NegBin = function(parameter, N, E, Ts){
      alpha = parameter[1]
      beta = parameter[2]

      max = sum(Ts *(dnbinom(N, size = alpha, prob = beta/(E + beta), log = TRUE)))

      return(-max)
    }

    theta_EM[i + 1, 1:2] = optim(c(theta_EM$alpha1[i], theta_EM$beta1[i]), NegBin, method = "BFGS", N = N, E = E, Ts = T_ij)$par
    theta_EM[i + 1, 3:4] = optim(c(theta_EM$alpha2[i], theta_EM$beta2[i]), NegBin, method = "BFGS", N = N, E = E, Ts = 1 - T_ij)$par

  }

  if (Loglik == TRUE){

    llh1 = rep(NA, nrow(theta_EM))

    for (i in c(1 : nrow(theta_EM))){

      probs1 <- theta_EM$beta1[i]/(Eij$baseline + theta_EM$beta1[i])
      probs2 <- theta_EM$beta2[i]/(Eij$baseline + theta_EM$beta2[i])

      ## Use LogSumExp trick
      D1.vec <- log(theta_EM$pi[i]) + dnbinom(Nij$frequency, size = theta_EM$alpha1[i], prob = probs1, log = TRUE)
      D2.vec <- log(1 - theta_EM$pi[i]) + dnbinom(Nij$frequency, size = theta_EM$alpha2[i], prob = probs2, log = TRUE)

      log.dens <- rep(0, length(D1.vec))
      log.dens[D1.vec < D2.vec] <- D2.vec[D1.vec < D2.vec] + log(1 + exp(D1.vec[D1.vec < D2.vec] - D2.vec[D1.vec < D2.vec]))
      log.dens[D1.vec >= D2.vec] <- D1.vec[D1.vec >= D2.vec] + log(1 + exp(D2.vec[D1.vec >= D2.vec] - D1.vec[D1.vec >= D2.vec]))

      llh1[i] = sum(log.dens)

    }
  } else {
    llh1 = NA
  }

  result = list("theta_EM" = theta_EM, "loglik" = llh1)
  return(result)

}






