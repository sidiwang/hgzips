#' Zero-Inflated Negative Binomial - Two Gamma Component
#'
#' This \code{ZINB_two_gamma} function finds hyperparameter estimates by implementing the Expectation-Maximization (EM) algorithm and zero-inflated negative binomial model with two gamma components. \code{nlminb} function is used to maximize the loglikelihood function.
#' @name ZINB_two_gamma
#' @aliases ZINB_two_gamma
#' @import stats
#'
#' @param alpha1 initial shape parameter value of the first gamma distribution for implementing the EM algprithm
#' @param beta1 initial rate parameter value of the first gamma distribution for implementing the EM algprithm
#' @param alpha2 initial shape parameter value of the second gamma distribution for implementing the EM algprithm
#' @param beta2 initial rate parameter value of the second gamma distribution for implementing the EM algprithm
#' @param pi initial mixing proportion guess of the two gamma distributions
#' @param omega initial weight for observing a true zero (according to zero-inflated poission distribution)
#' @param N vector of Nij values
#' @param E vector of Eij values
#' @param weight set weight = rep(1, length(N)) if N and E are not squashed data, or input the weight vector corresponding to the squashed Nij vector.
#' @param iteration number of EM algorithm iterations to run
#' @param zeroes A logical scalar specifying if zero counts should be included.
#' @param N_star the minimum Nij count size to be used for hyperparameter estimation. If zeroes are included in Nij vector, please set N_star = NULL
#' @param Loglik whether to return the loglikelihood of each iteration or not (TRUE or FALSE)
#'
#' @return a list of estimated parameters and their corresponding loglikelihood
#' if Nij = 0 is included in the input dataset
#' \itemize{
#' \item \code{theta_EM} Estimate of hyperparameters for each EM iteration
#' \item \code{llh} logliklihood for each EM iteration (optional)
#' }
#' if the minimum Nij count size included in the input dataset is not 0"
#' \itemize{
#' \item \code{theta_EM} Estimate of hyperparameters for each EM iteration
#' }
#'
#' @seealso openEBGM, nlminb
#' @examples
#' Nij = rnbinom(100, size = 2, prob = 0.3)
#' Eij = runif(100, 0, 2)
#' par_estimated = ZINB_two_gamma(0.2, 0.1, 2, 4, 0.33, 0.3, N = squashed$N, E = squashed$E, weight = squashed$weight, 10, Loglik = TRUE, zeroes = FALSE, N_star = 1)
#' par_estimated$theta_EM
#'
#' @export
#'

# input N = Nij$frequency, E = Eij$baseline


ZINB_two_gamma = function(alpha1, beta1, alpha2, beta2, pi, omega, N, E, weight, iteration, Loglik = FALSE, zeroes = FALSE, N_star = 1){

  squashedData = as.data.frame(cbind(N, E, weight))

  if (!is.null(N_star)){
    squashed = subset(squashed, N >= N_star)
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

  T_ij = rep(NA, nrow(squashedData))
  cluster = rep(NA, nrow(squashedData))
  u_ij = rep(NA, nrow(squashedData))


  # count = 0
  zero = which(N == 0)
  u_ij = ifelse(N == 0, 1, 0)

  # the E-M iterations
  for (i in c(1:N.EM)) {

    print(i)

    # Expectation step:
    if (zeroes == TRUE){
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

      # theta_EM[i + 1, 1:6] = optim(c(theta_EM$alpha1[i], theta_EM$beta1[i], theta_EM$alpha2[i], theta_EM$beta2[i], theta_EM$pi[i], theta_EM$omega[i]), ZNegBin, method = "BFGS", u = u_ij, N = N, E = E, weight = weight, Ts = T_ij)$par

      theta_EM[i + 1, 1:6] = stats::nlminb(c(theta_EM$alpha1[i], theta_EM$beta1[i], theta_EM$alpha2[i], theta_EM$beta2[i], theta_EM$pi[i], theta_EM$omega[i]), ZNegBin, u = u_ij, N = N, E = E, weight = weight, Ts = T_ij)$par


    }else{

      # Expectation step:
      probs1 <- theta_EM$beta1[i]/(squashed$E + theta_EM$beta1[i])
      probs2 <- theta_EM$beta2[i]/(squashed$E + theta_EM$beta2[i])
      T_ij_1 = dnbinom(squashed$N, size = theta_EM$alpha1[i], prob = probs1, log = TRUE) - log1p(-pnbinom((N_star - 1), size = theta_EM$alpha1[i], prob = probs1))
      T_ij_2 = dnbinom(squashed$N, size = theta_EM$alpha2[i], prob = probs2, log = TRUE) - log1p(-pnbinom((N_star - 1), size = theta_EM$alpha2[i], prob = probs2))

      logBF <- T_ij_2 - T_ij_1 + log(1 - theta_EM$pi[i]) - log(theta_EM$pi[i])
      T_ij[logBF < 0] <- exp(-log1p(exp(logBF[logBF < 0])))
      T_ij[logBF >= 0] <- exp(-logBF[logBF >= 0] - log1p(exp(-logBF[logBF>=0])))

      # Maximization step:

      theta_EM$pi[i + 1] = sum(T_ij*squashed$weight)/sum(squashed$weight)

      NegBin = function(parameter, N, E, W, Ts, zeroes){
        alpha = parameter[1]
        beta = parameter[2]

        max = sum(Ts*W*(dnbinom(N, size = alpha, prob = beta/(E + beta), log = TRUE) - log1p(-pnbinom((N_star - 1), size = alpha, prob = beta/(E + beta)))))

        return(-max)
      }

      theta_EM[i + 1, 1:2] = stats::nlminb(c(theta_EM$alpha1[i], theta_EM$beta1[i]), NegBin, N = squashed$N, E = squashed$E, W = squashed$weight, Ts = T_ij, zeroes = zeroes)$par
      theta_EM[i + 1, 3:4] = stats::nlminb(c(theta_EM$alpha2[i], theta_EM$beta2[i]), NegBin, N = squashed$N, E = squashed$E, W = squashed$weight, Ts = 1 - T_ij, zeroes = zeroes)$par

      N_ij_I = ifelse(N_ij == 0, 1, 0)
      v_.j = colSums(N_ij_I)
      b1 = theta_EM$beta1[i + 1]
      a1 = theta_EM$alpha1[i + 1]
      b2 = theta_EM$beta2[i + 1]
      a2 = theta_EM$alpha2[i + 1]
      pi = theta_EM$pi[i + 1]
      eta_j = rep(NA, ncol(N_ij))
      for (j in 1:ncol(N_ij)){
        eta_j[j] = sum(pi*(b1/((E_ij[, j] + b1))^a1) + (1 - pi)*(b2/((E_ij[, j] + b2)^a2)))
      }
      mean_eta_j = mean(eta_j)
      findroot = function(omega){
        return(nrow(N_ij)*omega + mean_eta_j - omega*mean_eta_j - mean(v_.j))
      }
      theta_EM$omega[i + 1] = uniroot(findroot, c(-1E6, 1E6), tol = 0.0001)$root
    }


  }



  if (Loglik == TRUE){

    llh1 = rep(NA, nrow(theta_EM))
    log.dens = rep(0, length(N))

    for (i in 1:nrow(theta_EM)){

      if (zeroes == TRUE){
        log.dens = ifelse(N == 0, log_sum_exp(log(theta_EM$omega[i]), (log(1 - theta_EM$omega[i]) + log_sum_exp((log(theta_EM$pi[i]) + dnbinom(0, size = theta_EM$alpha1[i], prob = theta_EM$beta1[i]/(E + theta_EM$beta1[i]), log = TRUE)), (log(1 - theta_EM$pi[i]) + dnbinom(0, size = theta_EM$alpha2[i], prob = theta_EM$beta2[i]/(E + theta_EM$beta2[i]), log = TRUE))))), (log(1 - theta_EM$omega[i]) + log_sum_exp((log(theta_EM$pi[i]) + dnbinom(N, size = theta_EM$alpha1[i], prob = theta_EM$beta1[i]/(E + theta_EM$beta1[i]), log = TRUE)), (log(1 - theta_EM$pi[i]) + dnbinom(N, size = theta_EM$alpha2[i], prob = theta_EM$beta2[i]/(E + theta_EM$beta2[i]), log = TRUE)))))
      }

      llh1[i] = sum(log.dens*weight)

    }
  } else {
    llh1 = NULL
  }

  result = list("theta_EM" = theta_EM, "Loglik" = llh1)
  return(result)

}





