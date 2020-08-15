#' HGZIPS - Data Squashing HZINB (not assuming independence)
#'
#' This Data Squashing HZINB function.........
#' @name gof
#' @import stats
#'
#' @param N 1
#' @param E 2
#' @param N_ij 3
#' @param E_ij 4
#' @param alpha1 5
#' @param alpha2 6
#' @param beta1 1
#' @param beta2 2
#' @param pi 3
#' @param omega 2
#' @param alpha 3
#' @param beta 4
#' @param posterior_a_j 2
#' @param posterior_b_j 2
#' @param posterior_omega_j 23
#' @param fitted_prob 3
#' @param observed_freq 4
#' @return observed frequencies, estimated frequencies and goodness of fit

#' @seealso
#'
##############
##########################################################
## goodness of fit
##########################################################

# input N_ij format
# output observed frequencies
#' @rdname gof
#' @return observed frequencies
#' observed_frequencies
#' @export
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

#' @rdname gof
#' @return ZINB fitted values
#' fitted_ZINB
#' @export
# input N = N_ij, E = E_ij
fitted_ZINB = function(alpha, beta, omega, N, E){

  N = as.matrix(N)
  E = as.matrix(E)

  fitted_prob_ZINB = ifelse(N == 0, omega + (1 - omega)*dnbinom(0, size = alpha, prob = beta/(E + beta)), (1 - omega)*dnbinom(N, size = alpha, prob = beta/(E + beta)))
  return(fitted_prob_ZINB)

}

#' @rdname gof
#' @return fitted MGPS
#' fitted MGPS
#' @export
# input N = N_ij, E = E_ij
fitted_MGPS = function(alpha1, alpha2, beta1, beta2, pi, N, E){

  N = as.matrix(N)
  E = as.matrix(E)

  fitted_prob_MGPS = pi*dnbinom(N, size = alpha1, prob = beta1/(E + beta1)) + (1 - pi)*dnbinom(N, size = alpha2, prob = beta2/(E + beta2))
  return(fitted_prob_MGPS)

}


#' @rdname gof
#' @return fitted ZINB two gamma
#' fitted ZINB two gamma
#' @export
# input N = N_ij, E = E_ij
fitted_ZINB_two_gamma = function(alpha1, alpha2, beta1, beta2, pi, omega, N, E){

  N = as.matrix(N)
  E = as.matrix(E)

  fitted_prob_ZINB_two_gamma = ifelse(N == 0, omega + (1 - omega)*(pi*dnbinom(0, size = alpha1, prob = beta1/(E + beta1)) + (1 - pi)*dnbinom(0, size = alpha2, prob = beta2/(E + beta2))), (1 - omega)*(pi*dnbinom(N, size = alpha1, prob = beta1/(E + beta1)) + (1 - pi)*dnbinom(N, size = alpha2, prob = beta2/(E + beta2))))
  return(fitted_prob_ZINB_two_gamma)

}


#' @rdname gof
#' @return observed_freq_HZINB
#' observed_freq_HZINB
#' @export
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




#' @rdname gof
#' @return fitted_HZINB_one_gamma_independence
#' fitted_HZINB_one_gamma_independence
#' @export
fitted_HZINB_one_gamma_independence = function(posterior_a_j, posterior_b_j, posterior_omega_j, N_ij, E_ij){

  N_ij = as.matrix(N_ij)
  E_ij = as.matrix(E_ij)

  fitted_prob_HZINB_one_gamma_independene = as.data.frame(matrix(NA, nrow(N_ij), ncol(N_ij)))

  for (j in 1:length(posterior_a_j)){

    fitted_prob_HZINB_one_gamma_independene[, j] = ifelse(N_ij[, j] == 0, posterior_omega_j[j] + (1 - posterior_omega_j[j])*(dnbinom(0, size = posterior_a_j[j], prob = posterior_b_j[j]/(E_ij[, j] + posterior_b_j[j]))), (1 - posterior_omega_j[j])*(pi*dnbinom(N_ij[,j], size = posterior_a_j[j], prob = posterior_b_j[j]/(E_ij[, j] + posterior_b_j[j]))))
  }

  return(fitted_prob_HZINB_one_gamma_independene)

}

#' @rdname gof
#' @return gof
#' gof
#' @export
gof = function(fitted_prob, observed_freq){

  gof = sum((fitted_prob - observed_freq)^2)
  return(gof)

}



