#' Multi-Item Gamma-Poisson Shrinker (MGPS)
#'
#' This \code{mgpsEM} function finds hyperparameter estimates by implementing the Expectation-Maximization (EM) algorithm.\code{nlminb} function is used to maximize the loglikelihood function.
#'
#' @name mgpsEM
#' @aliases mgpsEM
#' @import stats
#'
#' @param alpha1 initial shape parameter guess of the first gamma distribution
#' @param beta1 initial rate parameter guess of the first gamma distribution
#' @param alpha2 initial shape parameter guess of the second gamma distribution
#' @param beta2 initial rate parameter guess of the second gamma distribution
#' @param pi initial mixing proportion guess of the two gamma distributions
#' @param N vector of Nij values
#' @param E vector of Eij values
#' @param weight set weight = rep(1, length(N)) if N and E are not squashed data, or input the weight vector corresponding to the squashed Nij vector.
#' @param iteration number of EM algorithm iterations to run
#' @param Loglik whether to return the loglikelihood of each iteration or not (TRUE or FALSE)
#' @param zeroes A logical scalar specifying if zero counts should be included.
#' @param N_star the minimum Nij count size to be used for hyperparameter estimation. If zeroes are included in Nij vector, please set N_star = NULL
#'
#'
#' @return a list including the following:
#' \itemize{
#' \item \code{theta_EM} Estimate of hyperparameters for each EM iteration
#' \item \code{llh} logliklihood for each EM iteration (optional)
#' }
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency Tables, With an Application to the FDA Spontaneous Reporting System." The American Statistician, 53(3), 177-190.
#' @seealso openEBGM, nlminb
#' @examples
#' Nij = rnbinom(100, size = 2, prob = 0.3)
#' Eij = runif(100, 0, 2)
#' par_estimated = mgpsEM(2, 4, 0.2, 0.1, 0.33, Nij, Eij, rep(1, length(Nij)), iteration = 10, Loglik = TRUE, zeros = TRUE, N_star = 1)
#' par_estimated$theta_EM
#' par_estimated$Loglik
#'
#'
#' @export

#'
#----------------------------#
#-- (1) use optim function --#
#----------------------------#

mgpsEM = function(alpha1, beta1, alpha2, beta2, pi, N, E, weight, iteration, Loglik = FALSE, zeroes = FALSE, N_star = 1){

  squashed = as.data.frame(cbind(N, E, weight))

  if (!is.null(N_star)){
    squashed = subset(squashed, N >= N_star)
  }

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

    if (zeroes == FALSE){
      T_ij_1 = dnbinom(squashed$N, size = theta_EM$alpha1[i], prob = probs1, log = TRUE) - log1p(-pnbinom((N_star - 1), size = theta_EM$alpha1[i], prob = probs1))
      T_ij_2 = dnbinom(squashed$N, size = theta_EM$alpha2[i], prob = probs2, log = TRUE) - log1p(-pnbinom((N_star - 1), size = theta_EM$alpha2[i], prob = probs2))
    } else {
      T_ij_1 = dnbinom(squashed$N, size = theta_EM$alpha1[i], prob = probs1, log = TRUE)
      T_ij_2 = dnbinom(squashed$N, size = theta_EM$alpha2[i], prob = probs2, log = TRUE)
    }


    logBF <- T_ij_2 - T_ij_1 + log(1 - theta_EM$pi[i]) - log(theta_EM$pi[i])
    T_ij[logBF < 0] <- exp(-log1p(exp(logBF[logBF < 0])))
    T_ij[logBF >= 0] <- exp(-logBF[logBF >= 0] - log1p(exp(-logBF[logBF>=0])))

    # Maximization step:

    theta_EM$pi[i + 1] = sum(T_ij*squashed$weight)/sum(squashed$weight)

    NegBin = function(parameter, N, E, W, Ts, zeroes){
      alpha = parameter[1]
      beta = parameter[2]

      if (zeroes == FALSE){
        max = sum(Ts*W*(dnbinom(N, size = alpha, prob = beta/(E + beta), log = TRUE) - log1p(-pnbinom((N_star - 1), size = alpha, prob = beta/(E + beta)))))
      }else{
        max = sum(Ts*W*(dnbinom(N, size = alpha, prob = beta/(E + beta), log = TRUE)))
      }

      return(-max)
    }

   # theta_EM[i + 1, 1:2] = optim(c(theta_EM$alpha1[i], theta_EM$beta1[i]), NegBin, method = "Nelder-Mead", N = squashed$N, E = squashed$E, W = squashed$weight, Ts = T_ij, zeroes = zeroes)$par
  #  theta_EM[i + 1, 3:4] = optim(c(theta_EM$alpha2[i], theta_EM$beta2[i]), NegBin, method = "Nelder-Mead", N = squashed$N, E = squashed$E, W = squashed$weight, Ts = 1 - T_ij, zeroes = zeroes)$par
    theta_EM[i + 1, 1:2] = stats::nlminb(c(theta_EM$alpha1[i], theta_EM$beta1[i]), NegBin, N = squashed$N, E = squashed$E, W = squashed$weight, Ts = T_ij, zeroes = zeroes)$par
    theta_EM[i + 1, 3:4] = stats::nlminb(c(theta_EM$alpha2[i], theta_EM$beta2[i]), NegBin, N = squashed$N, E = squashed$E, W = squashed$weight, Ts = 1 - T_ij, zeroes = zeroes)$par
  }


  if (Loglik == TRUE){
    Tij = list()
    llh1 = rep(NA, nrow(theta_EM))

    for (i in c(1:nrow(theta_EM))){
     probs1 <- theta_EM$beta1[i]/(squashed$E + theta_EM$beta1[i])
     probs2 <- theta_EM$beta2[i]/(squashed$E + theta_EM$beta2[i])

     ## Use LogSumExp trick
     if (zeroes == FALSE){
       D1.vec <- log(theta_EM$pi[i]) + dnbinom(squashed$N, size = theta_EM$alpha1[i], prob = probs1, log = TRUE) - log1p(-pnbinom((N_star - 1), size = theta_EM$alpha1[i], prob = probs1))
       D2.vec <- log(1 - theta_EM$pi[i]) + dnbinom(squashed$N, size = theta_EM$alpha2[i], prob = probs2, log = TRUE) - log1p(-pnbinom((N_star - 1), size = theta_EM$alpha2[i], prob = probs2))
     }else{
       D1.vec <- log(theta_EM$pi[i]) + dnbinom(squashed$N, size = theta_EM$alpha1[i], prob = probs1, log = TRUE)
       D2.vec <- log(1 - theta_EM$pi[i]) + dnbinom(squashed$N, size = theta_EM$alpha2[i], prob = probs2, log = TRUE)
     }

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



