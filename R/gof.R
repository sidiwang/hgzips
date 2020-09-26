#' goodness-of-fit | fitted probabilities | observed probabilities
#'
#' Functions to calculate the goodness-of-fit, fitted probabilities, or observed probabilities of different models
#'
#'
#' @name gof
#' @import stats
#'
#' @param N vector of Nij values
#' @param E vector of Eij values
#' @param N_ij matrix of N_ij values, i - AE, j - drugs
#' @param E_ij matrix of E_ij values, corresponding to N_ij
#' @param posterior_a_j vector of posterior mean a_j values for each drug (j)
#' @param posterior_b_j vector of posterior mean b_j values for each drug (j)
#' @param posterior_omega_j vector of posterior mean omega_j values for each drug (j)
#' @param fitted_prob a dataset of calculated fitted probability for each AE and drug pair
#' @param observed_freq a dataset of calculated observed probability for each AE and drug pair
#' @details
#' Functions calculate observed probability for each AE and drug pair:
#' \code{observed_freqencies} for MGPS, ZINB one gamma and ZINB two gamma models
#' \code{observed_freq_HZINB} for HZINB model
#'
#' Functions calculate fitted probability for each AE and drug pair:
#' \code{fitted_MGPS} for MGPS model
#' \code{fitted_ZINB} for ZINB one gamma model
#' \code{fitted_ZINB_two_gamma} for ZINB two gamma model
#' \code{fitted_HZINB_one_gamma_independence} for HZINB model
#'
#' Function calculate goodness of fit of a model:
#' \code{gof}
#'
#' @seealso
#'
##############
##########################################################
## goodness of fit
##########################################################

# input N_ij format
# output observed frequencies
#' @rdname gof
#' @return \code{observed_freqencies} functions returns the observed probability for MGPS model.
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
#' @return \code{fitted_ZINB} functions returns the fitted probability for ZINB one gamma model.
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
#' @return \code{fitted_MGPS} functions returns the fitted probability for MGPS model.
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
#' @return \code{fitted_ZINB_two_gamma} functions returns the fitted probability for ZINB two gamma model.
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
#' @return \code{observed_freq_HZINB} functions returns the observed probability for HZINB model.
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
#' @return \code{fitted_HZINB_one_gamma_independence} functions returns the fitted probability for HZINB model.
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
#' @return \code{gof} functions returns the goodness-of-fit of a model.
#' gof
#' @export
gof = function(fitted_prob, observed_freq){

  gof = sum((fitted_prob - observed_freq)^2)
  return(gof)

}



