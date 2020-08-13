##########################################################
## Posterior HZINB One Gamma
##########################################################

# N = N_ij, E = E_ij

posteriorHZINB = function(posterior_a_j, posterior_b_j, posterior_omega_j, N_ij, E_ij){
  
  post_HZINB = matrix(NA, nrow(N_ij), ncol(N_ij))
  
  for (j in 1 : ncol(N_ij)){
    post_HZINB[, j] =  (posterior_a_j[j] + N_ij[, j])/(posterior_b_j[j] + E_ij[, j])
    Qx = posterior_omega_j[j]/(posterior_omega_j[j] + (1 - posterior_omega_j[j])*dnbinom(0, size = posterior_a_j[j], prob = posterior_b_j[j]/E_ij[, j] + posterior_b_j[j]))
    post_HZINB[which(N_ij[, j] == 0, arr.ind = TRUE)] = (Qx*posterior_a_j[j]/posterior_b_j[j] + (1 - Qx)*posterior_a_j[j]/(posterior_b_j[j] + E_ij[, j]))[which(N_ij[, j] == 0, arr.ind = TRUE)]
    
  } 
  
  return(post_HZINB)
}



