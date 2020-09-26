#' Zero-Inflated Negative Binomial - One Gamma Component
#'
#' This \code{ZINB_one_gamma} function finds hyperparameter estimates by implementing the Expectation-Maximization (EM) algorithm and zero-inflated negative binomial model with one gamma component. \code{nlminb} function is used to maximize the loglikelihood function.
#' @name ZINB_one_gamma
#' @aliases ZINB_one_gamma
#' @import stats
#'
#' @param alpha initial shape parameter value of the gamma distribution for implementing the EM algprithm
#' @param beta initial rate parameter value of the gamma distribution for implementing the EM algprithm
#' @param omega initial weight for observing a true zero (according to zero-inflated poission distribution)
#' @param N vector of Nij values
#' @param E vector of Eij values
#' @param weight set weight = rep(1, length(N)) if N and E are not squashed data, or input the weight vector corresponding to the squashed Nij vector.
#' @param iteration number of EM algorithm iterations to run
#' @param Loglik whether to return the loglikelihood of each iteration or not (TRUE or FALSE)
#' @param zeroes A logical scalar specifying if zero counts should be included.
#' @param N_star the minimum Nij count size to be used for hyperparameter estimation. If zeroes are included in Nij vector, please set N_star = NULL
#'
#' @return a list including the following:
#' if Nij = 0 is included in the input dataset
#' \itemize{
#' \item \code{theta_EM} Estimate of hyperparameters for each EM iteration
#' \item \code{llh} logliklihood for each EM iteration (optional)
#' }
#' if the minimum Nij count size included in the input dataset is not 0"
#' \itemize{
#' \item \code{alpha} Estimated alpha value
#' \item \code{beta} Estimated beta value
#' \item \code{omega} is not exported because the value of \code{omega} doesn't affect the estimation of \code{lambda}
#' }
#'
#' @seealso openEBGM, nlminb
#' @examples
#' Nij = rnbinom(100, size = 2, prob = 0.3)
#' Eij = runif(100, 0, 2)
#' par_estimated = ZINB_one_gamma(0.2, 0.1, 0.33, Nij, Eij, rep(1, length(Nij)), iteration = 10, Loglik = TRUE, zeros = TRUE, N_star = 1)
#' par_estimated$theta_EM
#' par_estimated$Loglik
#'
#'
#' @export
#'
ZINB_one_gamma = function(alpha, beta, omega, N, E, weight, iteration, Loglik = FALSE, zeroes = FALSE, N_star = 1){


  N_E_weight = as.data.frame(cbind(N, E, weight))
  N_E_weight$logisticY = ifelse(N_E_weight$N == 0, 0, 1)

  if (!is.null(N_star)){
    N_E_weight = subset(N_E_weight, N >= N_star)
  }

  if (zeroes == TRUE){

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

      ZNegBin = function(parameter, N, E, weight, Ts){
        alpha = parameter[1]
        beta = parameter[2]

        max = sum((1 - Ts)*weight*dnbinom(N, size = alpha, prob = beta/(E + beta), log = TRUE))

        return(-max)
      }

      #theta_EM[i + 1, 1:2] = optim(c(theta_EM$alpha[i], theta_EM$beta[i]), ZNegBin, method = "BFGS", N = N_E_weight$N, E = N_E_weight$E, weight = N_E_weight$weight, Ts = T_ij)$par
      theta_EM[i + 1, 1:2] =  stats::nlminb(c(theta_EM$alpha[i], theta_EM$beta[i]), ZNegBin, N = N_E_weight$N, E = N_E_weight$E, weight = N_E_weight$weight, Ts = T_ij)$par

    }

    if (Loglik == TRUE){

      llh1 = rep(NA, nrow(theta_EM))
      log.dens = rep(0, length(N_E_weight$N))

      for (i in c(1: nrow(theta_EM))){

        log.dens = ifelse(N_E_weight$N == 0, log(theta_EM$omega[i] + exp(log(1 - theta_EM$omega[i]) + dnbinom(0, size = theta_EM$alpha[i], prob = theta_EM$beta[i]/(N_E_weight$E + theta_EM$beta[i]), log = TRUE))), log(1 - theta_EM$omega[i]) + dnbinom(N_E_weight$N, size = theta_EM$alpha[i], prob = theta_EM$beta[i]/(N_E_weight$E + theta_EM$beta[i]), log = TRUE))

        llh1[i] = sum(log.dens*N_E_weight$weight)
      }
    } else {
      llh1 = NULL
    }

  }else{

    ZNegBin = function(parameter, N, E, weight){
      alpha = parameter[1]
      beta = parameter[2]
      omega = parameter[3]

      max = sum(weight*(dnbinom(N, size = alpha, prob = beta/(E + beta), log = TRUE) - log1p(-(pnbinom((N_star - 1), size = alpha, prob = beta/(E + beta))))))

      return(-max)
    }

    par =  stats::nlminb(c(theta_EM$alpha[i], theta_EM$beta[i]), ZNegBin, N = N_E_weight$N, E = N_E_weight$E, weight = N_E_weight$weight)$par
    alpha = par[1]
    beta = par[2]

    #   N_ij_I = ifelse(N_ij == 0, 1, 0)
    #    v_.j = colSums(N_ij_I)
    #  b = theta_EM$beta[i + 1]
    #  a = theta_EM$alpha[i + 1]
    #    eta_j = rep(NA, ncol(N_ij))
    #   for (j in 1:ncol(N_ij)){
    #    eta_j[j] = sum((b/(E_ij[, j] + b))^a)
    #  }
    #  mean_eta_j = mean(eta_j)
    #    findroot = function(omega){
    #     return(nrow(N_ij)*omega + mean_eta_j - omega*mean_eta_j - mean(v_.j))
    #  }
    #  theta_EM$omega[i + 1] = uniroot(findroot, c(-1E6, 1E6), tol = 0.0001)$root
  }

  if (zeroes == TRUE){

    result = list("theta_EM" = theta_EM, "Loglik" = llh1)

  }else{
    result = list("alpha" = alpha, "beta" = beta)
  }
  return(result)

}




