#' HGZIPS - MGPS (optim)
#'
#' This MGPS function.........
#' @name MGPS_optim
#' @aliases MGPS_optim
#' @import stats
#'
#' @param alpha1 initial shape parameter value of the first gamma distribution for implementing the EM algprithm
#' @param beta1 initial rate parameter value of the first gamma distribution for implementing the EM algprithm
#' @param alpha2 initial shape parameter value of the second gamma distribution for implementing the EM algprithm
#' @param beta2 initial rate parameter value of the second gamma distribution for implementing the EM algprithm
#' @param pi initial xxxxxxx?
#' @param N vector of Nij values
#' @param E vector of Eij values
#' @param weight set weight = rep(1, length(N)) if N and E are not squashed data, or input the weight vector corresponding to the squashed Nij vector.
#' @param iteration number of EM algorithm iterations to run
#' @param Loglik whether to return the loglikelihood of each iteration or not (TRUE or FALSE)
#' @return a list of estimated parameters and their corresponding loglikelihood
#' @export
#' @seealso
#'
#----------------------------#
#-- (1) use optim function --#
#----------------------------#

MGPS_optim = function(alpha1, beta1, alpha2, beta2, pi, N, E, weight, iteration, Loglik){

  squashed = as.data.frame(cbind(N, E, weight))

  # initialization
  N.EM <- iteration  # number of E-M iterations
  theta_EM = matrix(rep(0, 5 * (N.EM + 1)), nrow = N.EM + 1, ncol = 5)  # matrix for holding theta in the E-M iterations
  theta_EM = as.data.frame(theta_EM)
  theta_EM[1, ] <- c(alpha1, beta1, alpha2, beta2, pi)  # initial value for parameters
  colnames(theta_EM) = c("alpha1", "beta1", "alpha2", "beta2", "pi")

  T_ij = rep(NA, nrow(squashed))
  cluster = rep(NA, nrow(squashed))


  # the E-M iterations
  for (i in c(1:N.EM)) {
    cat("\nIteration", i, "\n")
    print(theta_EM[i, ])

    # Expectation step:
    probs1 <- theta_EM$beta1[i]/(squashed$E + theta_EM$beta1[i])
    probs2 <- theta_EM$beta2[i]/(squashed$E + theta_EM$beta2[i])
    T_ij_1 = dnbinom(squashed$N, size = theta_EM$alpha1[i], prob = probs1, log = TRUE)
    T_ij_2 = dnbinom(squashed$N, size = theta_EM$alpha2[i], prob = probs2, log = TRUE)

    logBF <- T_ij_2 - T_ij_1 + log(1 - theta_EM$pi[i]) - log(theta_EM$pi[i])
    T_ij[logBF < 0] <- exp(-log1p(exp(logBF[logBF < 0])))
    T_ij[logBF >= 0] <- exp(-logBF[logBF >= 0] - log1p(exp(-logBF[logBF>=0])))

    # Maximization step:

    theta_EM$pi[i + 1] = sum(T_ij*squashed$weight)/sum(squashed$weight)

    NegBin = function(parameter, N, E, W, Ts){
      alpha = parameter[1]
      beta = parameter[2]

      max = sum(Ts*W*(dnbinom(N, size = alpha, prob = beta/(E + beta), log = TRUE)))

      return(-max)
    }

    theta_EM[i + 1, 1:2] = optim(c(theta_EM$alpha1[i], theta_EM$beta1[i]), NegBin, method = "BFGS", N = squashed$N, E = squashed$E, W = squashed$weight, Ts = T_ij)$par
    theta_EM[i + 1, 3:4] = optim(c(theta_EM$alpha2[i], theta_EM$beta2[i]), NegBin, method = "BFGS", N = squashed$N, E = squashed$E, W = squashed$weight, Ts = 1 - T_ij)$par

  }


  if (Loglik == TRUE){
    Tij = list()
    llh1 = rep(NA, nrow(theta_EM))

    for (i in c(1:nrow(theta_EM))){
     probs1 <- theta_EM$beta1[i]/(squashed$E + theta_EM$beta1[i])
     probs2 <- theta_EM$beta2[i]/(squashed$E + theta_EM$beta2[i])

     ## Use LogSumExp trick
     D1.vec <- log(theta_EM$pi[i]) + dnbinom(squashed$N, size = theta_EM$alpha1[i], prob = probs1, log = TRUE)
     D2.vec <- log(1 - theta_EM$pi[i]) + dnbinom(squashed$N, size = theta_EM$alpha2[i], prob = probs2, log = TRUE)
     log.dens <- rep(0, length(D1.vec))
     log.dens[D1.vec < D2.vec] <- D2.vec[D1.vec < D2.vec] + log(1 + exp(D1.vec[D1.vec < D2.vec] - D2.vec[D1.vec < D2.vec]))
     log.dens[D1.vec >= D2.vec] <- D1.vec[D1.vec >= D2.vec] + log(1 + exp(D2.vec[D1.vec >= D2.vec] - D1.vec[D1.vec >= D2.vec]))

     llh1[i] =  sum(log.dens*squashed$weight)
     }
    } else {
    llh1 = NULL
  }

  result = list("theta_EM" = theta_EM, "Loglik" = llh1)
  return(result)
}



