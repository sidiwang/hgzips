###########################################################
## Posterior expectation of lambda|N (MGPS)
###########################################################

# input, N = Nij$frequency, E = Eij$baseline

post_mean_lambda_ZINB = function(alpha, beta, N, E){
  post_mean_lambda = (alpha + Nij$frequency)/(beta + Eij$baseline)
  return(post_mean_lambda)
}

post_mean_loglambda_ZINB = function(alpha, beta, N, E){
  post_mean_loglambda = digamma(alpha + Nij$frequency) - log(beta + Eij$baseline)
  return(post_mean_loglambda)
}