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
#'
##########################################################
## goodness of fit
##########################################################

# input N_ij format
# output observed frequencies

observed_freqencies = function(N){

  N = as.matrix(N)

  observed_freq = as.data.frame(matrix(NA, nrow(N), ncol(N)))
  N_ij = as.data.frame(N)

  aa = unique(as.vector(as.matrix(N_ij)))

  for (i in aa){
    observed_freq[which(N_ij == i, arr.ind = TRUE)] = length(which(N_ij == i))/length(as.vector(as.matrix(N_ij)))
  }

  return(observed_freq)

}


# input N = N_ij, E = E_ij
fitted_ZINB = function(alpha, beta, omega, N, E){

  N = as.matrix(N)
  E = as.matrix(E)

  fitted_prob_ZINB = ifelse(N == 0, omega + (1 - omega)*dnbinom(0, size = alpha, prob = beta/(E + beta)), (1 - omega)*dnbinom(N, size = alpha, prob = beta/(E + beta)))
  return(fitted_prob_ZINB)

}


# input N = N_ij, E = E_ij
fitted_MGPS = function(alpha1, alpha2, beta1, beta2, pi, N, E){

  N = as.matrix(N)
  E = as.matrix(E)

  fitted_prob_MGPS = pi*dnbinom(N, size = alpha1, prob = beta1/(E + beta1)) + (1 - pi)*dnbinom(N, size = alpha2, prob = beta2/(E + beta2))
  return(fitted_prob_MGPS)

}

# input N = N_ij, E = E_ij
fitted_ZINB_two_gamma = function(alpha1, alpha2, beta1, beta2, pi, omega, N, E){

  N = as.matrix(N)
  E = as.matrix(E)

  fitted_prob_ZINB_two_gamma = ifelse(N == 0, omega + (1 - omega)*(pi*dnbinom(0, size = alpha1, prob = beta1/(E + beta1)) + (1 - pi)*dnbinom(0, size = alpha2, prob = beta2/(E + beta2))), (1 - omega)*(pi*dnbinom(N, size = alpha1, prob = beta1/(E + beta1)) + (1 - pi)*dnbinom(N, size = alpha2, prob = beta2/(E + beta2))))
  return(fitted_prob_ZINB_two_gamma)

}


observed_freq_HZINB = function(N_ij){

  observed_freq = as.data.frame(matrix(NA, nrow(N_ij), ncol(N_ij)))

  N_ij = as.data.frame(N_ij)

  for (i in 1:nrow(N_ij)){
    for (j in 1:ncol(N_ij)){
      observed_freq[i,j] = length(which(N_ij[,j] == N_ij[i,j]))/length(N_ij[,j])
    }
  }

  return(observed_freq)

}





fitted_HZINB_one_gamma_independence = function(posterior_a_j, posterior_b_j, posterior_omega_j, N_ij, E_ij){

  N_ij = as.matrix(N_ij)
  E_ij = as.matrix(E_ij)

  fitted_prob_HZINB_one_gamma_independene = as.data.frame(matrix(NA, nrow(N_ij), ncol(N_ij)))

  for (j in 1:length(posterior_a_j)){

    fitted_prob_HZINB_one_gamma_independene[, j] = ifelse(N_ij[, j] == 0, posterior_omega_j[j] + (1 - posterior_omega_j[j])*(dnbinom(0, size = posterior_a_j[j], prob = posterior_b_j[j]/(E_ij[, j] + posterior_b_j[j]))), (1 - posterior_omega_j[j])*(pi*dnbinom(N_ij[,j], size = posterior_a_j[j], prob = posterior_b_j[j]/(E_ij[, j] + posterior_b_j[j]))))
  }

  return(fitted_prob_HZINB_one_gamma_independene)

}

gof = function(fitted_prob, observed_freq){

  gof = sum((fitted_prob - observed_freq)^2)
  return(gof)

}



