##########################################################
## HZINB: check the range of parameters
##########################################################

posterior_abomega = function(grid, pi_klh_final_a_j, pi_klh_final_b_j, pi_klh_final_omega_j, N_ij, E_ij){
  
  a_j = grid$a_j
  b_j = grid$b_j
  omega_j = grid$omega_j
  
  K = length(grid$a_j)
  L = length(grid$b_j)
  H = length(grid$omega_j)
  
  all_combinations = as.data.frame(matrix(NA, K*L*H, 3))
  colnames(all_combinations) = c("a_j", "b_j", "omega_j")
  all_combinations$a_j = rep(grid$a_j, 100)
  
  for (i in c(1:L)){
    all_combinations$b_j[((i - 1)*100 + 1):(i*100)] = rep(grid$b_j[i], 100)
  }
  
  for (i in c(1:H)){
    row_num = c(((i - 1)*10 + 1):(i*10))
    all_combinations$omega_j[row_num] = rep(grid$omega_j[i], 10)
  }
  
  all_combinations$omega_j = rep(all_combinations$omega_j[1:100], 10)
  
  joint_probs = as.data.frame(matrix(NA, nrow(all_combinations), ncol(N_ij)))
  
  for (j in 1:ncol(N_ij)){
    for (m in 1:nrow(all_combinations)){
      joint_probs[m,j] = sum(dzinbinom(N_ij[,j], mu = (E_ij[,j]/all_combinations$b_j[m])*all_combinations$a_j[m], size = all_combinations$a_j[m], pi = all_combinations$omega_j[m], log = TRUE))
    }
  }
  
  pi_klh_all_combinations = as.data.frame(matrix(NA, K*L*H, 3))
  colnames(pi_klh_all_combinations) = c("K", "L", "H")
  pi_klh_all_combinations$K = rep(pi_klh_final_a_j[1,], 100)
  
  for (i in c(1:L)){
    pi_klh_all_combinations$L[((i - 1)*100 + 1):(i*100)] = rep(pi_klh_final_b_j[1,][i], 100)
  }
  
  for (i in c(1:H)){
    row_num = c(((i - 1)*10 + 1):(i*10))
    pi_klh_all_combinations$H[row_num] = rep(pi_klh_final_omega_j[1,][i], 10)
  }
  
  pi_klh_all_combinations$H = rep(pi_klh_all_combinations$H[1:100], 10)
  
  pi_klh_all_combinations$prod = unlist(pi_klh_all_combinations$K) * unlist(pi_klh_all_combinations$L) * unlist(pi_klh_all_combinations$H)
  
  LSE_R <- function(vec){ 
    n.vec <- length(vec)
    vec <- sort(vec, decreasing = TRUE)
    Lk <- vec[1]
    for (k in 1:(n.vec-1)) {
      Lk <- max(vec[k+1], Lk) + log1p(exp(-abs(vec[k+1] - Lk))) 
    }
    return(Lk)
  }
  
  expectation_a_j = rep(NA, ncol(N_ij))
  expectation_b_j = rep(NA, ncol(N_ij))
  expectation_omega_j = rep(NA, ncol(N_ij))
  
  for (j in 1 : ncol(N_ij)){
    
    print(j)
    
    denominator_a_j = LSE_R(joint_probs[, j] + log(pi_klh_all_combinations$prod))
    numerator_a_j = LSE_R(joint_probs[, j] + log(all_combinations$a_j) + log(pi_klh_all_combinations$prod))
    expectation_a_j[j] = exp(numerator_a_j - denominator_a_j)
    
    denominator_b_j = LSE_R(joint_probs[, j] + log(pi_klh_all_combinations$prod))
    numerator_b_j = LSE_R(joint_probs[, j] + log(all_combinations$b_j) + log(pi_klh_all_combinations$prod))
    expectation_b_j[j] = exp(numerator_b_j - denominator_b_j)
    
    denominator_omega_j = LSE_R(joint_probs[, j] + log(pi_klh_all_combinations$prod))
    numerator_omega_j = LSE_R(joint_probs[, j] + log(all_combinations$omega_j) + log(pi_klh_all_combinations$prod))
    expectation_omega_j[j] = exp(numerator_omega_j - denominator_omega_j)
    
  }
  
  result = list("posterior_a_j" = expectation_a_j, "posterior_b_j" = expectation_b_j, "posterior_omega_j" = expectation_omega_j)
  return(result)
  
}
