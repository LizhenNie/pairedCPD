# this function is for the case where we use Algorithm I, and the two naive baselines use RBF kernel
all_stats = function(x, y, hx, n0, n1, h_original = 0.1, h_concate = 0.1){
  
  # x is n * p matrix
  # y is n * 1 matrix
  
  K = exp( - as.matrix(dist(x)^2) / hx)
  n = length(y)
  K_row_cusum = t(apply(K, 1, cumsum)) # rowwise cumulative sum
  Ky = K %*% diag(c(y)) # scale each column of K by y[i]
  Ky_row_cusum = t(apply(Ky, 1, cumsum))
  
  # ------------------- two one-sided estimators ------------------
  delta = rep(0, n)
  mdelta = rep(0, n)
  
  for (i in n0:n1){
    
    # improved estimator (at every point x0): random design
    tmp_left = Ky_row_cusum[, i] / K_row_cusum[, i]
    tmp_right = (Ky_row_cusum[, n] - Ky_row_cusum[, i]) / (K_row_cusum[, n] - K_row_cusum[, i])
    
    mdelta[i] = mean((tmp_left - tmp_right)^2, na.rm = T) # all x0's
    delta[i] = tmp_left[i] - tmp_right[i] # random x0
  }
  
  # random design (asymptotics rescaling)
  myseq = 1 : n 
  
  # at every point x0
  all_seq = mdelta * ((myseq - 1) * (n - myseq)) / n ^ 2
  
  # -------------- one function fits all --------------
  fitted = Ky_row_cusum[, n] / K_row_cusum[, n]
  delta_10 = y - fitted # get residuals
  
  tmp = our_method(x = delta_10, n0 = n0, n1 = n1)
  tmp_one_all = tmp$Z1
  loc_10 = as.numeric(tmp$location1)
  
  # -------------- one function fits all (original data) ----------------
  K = exp( - as.matrix(dist(x)^2) / h_original)
  tmp = our_method_K(K = K, n0 = n0, n1 = n1)
  tmp_one_all_original = tmp$Z1
  loc_10_ori = as.numeric(tmp$location1)
  
  # -------------- one function fits all (concatenated data) ----------------
  x = as.matrix(x)
  K = exp( - as.matrix(dist(cbind(y, x))^2) / h_concate)
  tmp = our_method_K(K = K, n0 = n0, n1 = n1)
  tmp_one_all_concate = tmp$Z1
  loc_concate = as.numeric(tmp$location1)
  
  return(list(two_sided_all = max(all_seq, na.rm = T),
              one_fits_all = tmp_one_all, 
              one_fits_all_original = tmp_one_all_original,
              one_concate = tmp_one_all_concate,
              two_sided_all_loc = which.max(all_seq),
              one_fits_all_loc = loc_10, 
              one_fits_all_original_loc = loc_10_ori,
              one_concate_loc = loc_concate
  ))
  
}

