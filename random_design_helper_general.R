# this function is for the general problem (conditional distribution changes)
# and the two naive baselines use RBF kernel
all_stats_general = function(Kx_res, Ky_res, Kx_two, Ky_two, K_concate, Ky_only, n0, n1){
  
  # x is n * p matrix
  # y is n * d matrix
  
  library(tictoc)
  
  K = Kx_two
  n = nrow(K)
  K_row_cusum = t(apply(K, 1, cumsum)) # rowwise cumulative sum
  KKT = K %*% K
  # Ky = y %*% t(y)
  K_twosides = KKT * Ky_two / n
  
  delta = rep(0, n)
  mdelta = rep(0, n)
  residual_delta = rep(0, n)
  
  # ------------------- two-sided estimators ------------------
  tic()
  for (i in n0:n1){
    
    # two-sided estimator (at every point x0)
    weight1 = 1 / K_row_cusum[1:i, i]
    weight2 = 1 / (K_row_cusum[(i + 1):n, n] - K_row_cusum[(i + 1):n, i])
    weights = c(weight1, weight2)
    reweightd_K_twosides = t(t(K_twosides * weights) * weights)
    mdelta[i] = sum(reweightd_K_twosides[1:i, 1:i]) + sum(reweightd_K_twosides[(i + 1):n, (i + 1):n]) - 2 * sum(reweightd_K_twosides[1:i, (i + 1):n]) # all x0's
    
  }
  toc()
  
  K = Kx_res
  n = nrow(K)
  K_row_cusum = t(apply(K, 1, cumsum)) # rowwise cumulative sum
  
  # ------------- residual-based estimator ---------------
  tic()
  for (i in n0:n1){
    
    weight1 = n / (n - i) * (K_row_cusum[1:i, n] - K_row_cusum[1:i, i]) / K_row_cusum[1:i, n]
    weight2 = n / i * K_row_cusum[(i + 1):n, i] / K_row_cusum[(i + 1):n, n]
    weights = c(weight1, weight2)
    
    reweighted_Ky = t(t(Ky_res * weights) * weights)
    residual_delta[i] = mean(reweighted_Ky[1:i, 1:i]) + mean(reweighted_Ky[(i + 1):n, (i + 1):n]) - 2 * mean(reweighted_Ky[1:i, (i + 1):n])
  }
  toc()
  
  # random design (asymptotics rescaling)
  myseq = 1 : n 
  
  # at every point x0
  all_seq = mdelta * ((myseq - 1) * (n - myseq)) / n ^ 2
  tmp_one_all = residual_delta * ((myseq - 1) * (n - myseq)) / n ^ 2
  
  # -------------- naive baseline I (original data) ----------------
  tic()
  tmp = our_method_K(K = Ky_only, n0 = n0, n1 = n1)
  tmp_one_all_original = tmp$Z1
  loc_10_ori = as.numeric(tmp$location1)
  toc()
  
  # -------------- naive baseline II (concatenated data) ----------------
  tic()
  tmp = our_method_K(K = K_concate, n0 = n0, n1 = n1)
  tmp_one_all_concate = tmp$Z1
  loc_concate = as.numeric(tmp$location1)
  toc()
  
  return(list(two_sided_all = max(all_seq, na.rm = T),
              one_fits_all = max(tmp_one_all, na.rm = T), 
              one_fits_all_original = tmp_one_all_original,
              one_concate = tmp_one_all_concate,
              two_sided_all_loc = which.max(all_seq),
              one_fits_all_loc = which.max(tmp_one_all), 
              one_fits_all_original_loc = loc_10_ori,
              one_concate_loc = loc_concate
  ))
  
}
