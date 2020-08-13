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

ProfileLogLik <- function(alpha, Tij, N, E) {

  N = as.matrix(N)
  E = as.matrix(E)

  fff <- function(x) {
    ## Think of x as beta that you want to find
    tau <- length(x)
    ans <- x
    for(k in 1:tau) {
      ans[k] <- (alpha/x[k])*sum(Tij) - sum((Tij*(N + alpha))/(x[k] + E))
    }
    return(ans)
  }
  beta.val <- uniroot(fff, interval = c(0, 1000), extendInt = "yes")$root

  probs <- beta.val/(E + beta.val)
  ans <- sum(Tij*dnbinom(N, size = alpha, prob = probs, log = TRUE))
  return(ans)
}


MGPSParUpdate <- function(alpha.vec, beta.vec, pi.vec, N, E) {

  N = as.matrix(N)
  E = as.matrix(E)

  I <- nrow(N)
  J <- ncol(N)

  LT1ij <- LT2ij <- post.probs <- rep(0, I)

  ff <- function(x, Tij, alpha) {
    tau <- length(x)
    ans <- x
    for(k in 1:tau) {
      ans[k] <- (alpha/x[k])*sum(Tij) - sum((Tij*(N + alpha))/(x[k] + E))
    }
    return(ans)
  }

  LT1ij <- dnbinom(N, size = alpha.vec[1], prob=beta.vec[1]/(E + beta.vec[1]), log=TRUE)
  LT2ij <- dnbinom(N, size = alpha.vec[2], prob=beta.vec[2]/(E + beta.vec[2]), log=TRUE)

  logBF <- LT2ij - LT1ij + log(pi.vec[2]) - log(pi.vec[1])
  post.probs[logBF < 0] <- exp(-log1p(exp(logBF[logBF < 0])))
  post.probs[logBF >= 0] <- exp(-logBF[logBF >= 0] - log1p(exp(-logBF[logBF>=0])))

  pi.vec <- c( mean(c(post.probs), na.rm = TRUE), 1 - mean(post.probs, na.rm = TRUE))

  alpha.vec[1] <- optimize(ProfileLogLik, c(0, 1000), maximum = TRUE, Tij = post.probs, N = N, E = E)$maximum
  beta.vec[1] <- uniroot(ff, interval=c(0, 1000), Tij = post.probs, alpha = alpha.vec[1], extendInt = "yes")$root

  alpha.vec[2] <- optimize(ProfileLogLik, c(0, 1000), maximum = TRUE, Tij = 1 - post.probs, N = N, E = E)$maximum
  beta.vec[2] <- uniroot(ff, interval = c(0, 1000), Tij = 1 - post.probs, alpha = alpha.vec[2], extendInt = "yes")$root

  par.list <- list("alpha1" = alpha.vec[1], "beta1" = beta.vec[1], "alpha2" = alpha.vec[2], "beta2" = beta.vec[2], "pi" = pi.vec[1])
  return(par.list)
}

LogLikMGPS <- function(theta, N, E) {
  N = as.matrix(N)
  E = as.matrix(E)
  I <- nrow(N)
  J <- ncol(N)
  D1 <- D2 <- matrix(0, nrow=I, ncol=J)
  for(k in 1:I) {
    D1[k,] <- dnbinom(N[k,], size = theta$alpha1, prob=theta$beta1/(E[k,] + theta$beta1), log=TRUE)
    D2[k,] <- dnbinom(N[k,], size = theta$alpha2, prob=theta$beta2/(E[k,] + theta$beta2), log=TRUE)
  }

  ## Use LogSumExp trick
  D1.vec <- log(theta$pi) + c(D1)
  D2.vec <- log(1 - theta$pi) + c(D2)
  log.dens <- rep(0, length(D1.vec))
  log.dens[D1.vec < D2.vec] <- D2.vec[D1.vec < D2.vec] + log(1 + exp(D1.vec[D1.vec < D2.vec] - D2.vec[D1.vec < D2.vec]))
  log.dens[D1.vec >= D2.vec] <- D1.vec[D1.vec >= D2.vec] + log(1 + exp(D2.vec[D1.vec >= D2.vec] - D1.vec[D1.vec >= D2.vec]))

  loglik <- sum(log.dens)
  return(loglik)
}


# input: initial alpha, beta, pi,  N = Nij$frequency, E = Eij$baseline, iterations
# output: estimated parameters of each iteration, and loglikelihood of each iteration

MGPS_uniroot_optimization = function(alpha.par, beta.par, pi.par, N, E, iterations, Loglik){

  niter <- iterations
  alpha.par = matrix(NA, 2, niter + 1)
  alpha.par[ ,1] = c(0.2, 2)
  beta.par = matrix(NA, 2, niter + 1)
  beta.par[, 1] = c(0.1, 4)
  pi.par = matrix(NA, 2, niter + 1)
  pi.par[, 1] = c(1/3, 2/3)
  theta0 <- list()
  theta0$alpha1 <- alpha.par[1,1]
  theta0$alpha2 <- alpha.par[2,1]
  theta0$beta1 <- beta.par[1,1]
  theta0$beta2 <- beta.par[2,1]
  theta0$pi <- pi.par[1,1]

  ell <- rep(0, niter + 1)
  ell[1] <- LogLikMGPS(theta0, N = N, E = E)

  for (i in 1:niter) {

    theta_EM = MGPSParUpdate(alpha.vec = alpha.par[, i], beta.vec = beta.par[, i], pi.vec = pi.par[, i], N = N, E = E)
    alpha.par[, i+1] = c(theta_EM$alpha1, theta_EM$alpha2)
    beta.par[, i+1] = c(theta_EM$beta1, theta_EM$beta2)
    pi.par[, i+1] = c(theta_EM$pi, 1 - theta_EM$pi)

    if (Loglik == TRUE){
      ell[i+1] <- LogLikMGPS(theta_EM, N = N, E = E)
    } else {
      ell = NA
    }
  }

  result = list("alpha" = alpha.par, "beta" = beta.par, "pi" = pi.par, "loglik" = ell)
  return(result)

}
