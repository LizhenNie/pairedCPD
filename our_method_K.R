our_method = function(x = x, n0 = NULL, n1 = NULL, correct = TRUE){
  
  library(RSpectra)
  library(tictoc)
  
  d <- as.matrix(dist(x)^2)
  n = nrow(d)
  for (i in 1:n){
    d[i,i] = 0
  }
  if(is.null(n0)){n0 = ceiling(0.05 * n)}
  if(is.null(n1)){n1 = floor(0.95 * n)}
  
  # quantities used in S1
  colsum = apply(d, 1, sum)
  totalsum = sum(colsum)
  s = sum((colsum * 2 / n - totalsum / n^2)^2) / 4 / n - 1 / 4 * totalsum^2 / n^4
  s = sqrt(s)
  dd = matrix(0, nrow = n + 1, ncol = n + 1)
  for (row in 2:(n + 1)){
    for (col in 2:(n + 1)){
      dd[row,col] = dd[row,col-1] + dd[row-1, col] - dd[row-1, col-1] + d[row-1, col-1]
    }
  }
  dd = dd[2:(n+1), 2:(n+1)]
  b1 = function(t){return(dd[t,t] / t^2)}
  b2 = function(t){ return( ( dd[n,n] - dd[t,n] - dd[n,t] + dd[t,t] ) / (n-t)^2 ) }
  a = function(t){ return( ( dd[n,t] - dd[t,t] ) / t / (n-t) ) }
  
  # calculating T1, T2 and T3
  bb1 = rep(0, n); bb2 = rep(0, n); aa = rep(0, n); 
  S1_f = rep(0, n); S2_f = rep(0, n); S3_f = rep(0, n); 
  correct_S2_f = rep(0, n)
  e2 = dd[n, n] / 2 / n / (n-1)
  
  for (t in n0:n1){
    bb1[t] = b1(t)
    bb2[t] = b2(t)
    aa[t] = a(t)
    if(correct){
      S1_f[t] = t * (n-t) / n * ( aa[t] - 0.5 * bb1[t] * t / (t-1) - 0.5 * bb2[t] * (n-t) / (n-t-1) )
    }else{
      S1_f[t] = t * (n-t) / n * ( aa[t] - 0.5 * bb1[t] - 0.5 * bb2[t] ) 
    }
    S2_f[t] = 1 / s / 2 * sqrt(t * (n-t) / n) * (bb1[t] - bb2[t])
    S3_f[t] = 1 / s^2 * n / t / (n-t) * S1_f[t]^2  + (S2_f[t] - 1 / sqrt(n) / s / sqrt( t/n * (1-t/n) ) * e2 * (t/n * (n-t)/(n-t-1) - t/(t-1 ) *(n-t)/n) )^2
    # S3 has to use same distance for S1 and S2, and we choose d = || ||^2 here, and by default S3 is using higher order corrections 
    
    if (correct){# using this correction has significant gains in high dimensions (both p-value calibration and localization)
      correct_S2_f[t] = abs(S2_f[t] - 1 / sqrt(n) / s / sqrt( t/n * (1-t/n) ) * e2 * (t/n*(n-t)/(n-t-1) - t/(t-1)*(n-t)/n) )
    }else{
      S2_f[t] = abs(S2_f[t])
    }
    
  }
  
  ####### calculating S1, S2, S3 ########
  Z1 = max(S1_f[n0:n1])
  predicted1 = which.max(S1_f[n0:n1]) + n0-1
  
  if(!correct){ # if does not do correction
    Z2 = max(S2_f[n0:n1])
    predicted2 = which.max(S2_f[n0:n1]) + n0-1
  }else{
    Z2 = max(correct_S2_f[n0:n1])
    predicted2 = which.max(correct_S2_f[n0:n1]) + n0-1
  }
  
  Z3 = max(S3_f[n0:n1])
  predicted3 = which.max(S3_f[n0:n1]) + n0-1
  predictedDubey = which.max(D_f[n0:n1]) + n0-1

  return(list(Z1 = Z1, location1 = predicted1,  
              Z2 = Z2, location2 = predicted2, 
              Z2 = Z3, location3 = predicted3))
  
}

our_method_K = function(K, n0 = NULL, n1 = NULL, correct = TRUE){
  
  # this one uses different distances for S1 and S2; for S3, we should use the same distance for performance guarantees
  
  library(RSpectra)
  library(tictoc)
  
  d = K
  
  n = nrow(d)
  if(is.null(n0)){n0 = ceiling(0.05 * n)}
  if(is.null(n1)){n1 = floor(0.95 * n)}
  
  # quantities used in S1
  colsum = apply(d, 1, sum)
  totalsum = sum(colsum)
  s = sum((colsum * 2 / n - totalsum / n^2)^2) / 4 / n - 1 / 4 * totalsum^2 / n^4
  s = sqrt(s)
  dd = matrix(0, nrow = n + 1, ncol = n + 1)
  for (row in 2:(n+1)){
    for (col in 2:(n+1)){
      dd[row,col] = dd[row, col-1] + dd[row-1, col] - dd[row-1, col-1] + d[row-1, col-1]
    }
  }
  dd = dd[2:(n+1), 2:(n+1)]
  b1 = function(t){return(dd[t,t] / t^2)}
  b2 = function(t){ return( ( dd[n,n] - dd[t,n] - dd[n,t] + dd[t,t] ) / (n-t)^2 ) }
  a = function(t){ return( ( dd[n,t] - dd[t,t] ) / t / (n-t) ) }
  
  # calculating T1, T2 and T3
  bb1 = rep(0, n); bb2 = rep(0, n); aa = rep(0, n)
  S1_f = rep(0, n); S2_f = rep(0, n); S3_f = rep(0, n); 
  correct_S2_f = rep(0, n)
  e2 = dd[n, n] / 2 / n / (n-1)
  
  for (t in n0:n1){
    
    bb1[t] = b1(t)
    bb2[t] = b2(t)
    aa[t] = a(t)
    
    if(correct){
      S1_f[t] = t * (n-t) / n * ( aa[t] - 0.5 * bb1[t] * t/(t-1) - 0.5 * bb2[t] * (n-t)/(n-t-1) )
    }else{
      S1_f[t] = t * (n-t) / n * ( aa[t] - 0.5 * bb1[t] - 0.5 * bb2[t] ) 
    }
    S2_f[t] = 1 / s / 2 * sqrt(t * (n-t) / n) * (bb1[t] - bb2[t])
    S3_f[t] = 1 / s^2 * n / t / (n-t) * S1_f[t]^2  + (S2_f[t] - 1 / sqrt(n) / s / sqrt( t/n * (1-t/n) ) * e2 * (t/n*(n-t)/(n-t-1) - t/(t-1)*(n-t)/n) )^2
    # S3 has to use same distance for S1 and S2, and we choose d = || ||^2 here, and by default S3 is using higher order corrections 
    
    if (correct){# using this correction has significant gains in high dimensions (both p-value calibration and localization)
      correct_S2_f[t] = abs(S2_f[t] - 1 / sqrt(n) / s / sqrt( t/n * (1-t/n) ) * e2 * (t/n*(n-t)/(n-t-1) - t/(t-1)*(n-t)/n) )
    }else{
      S2_f[t] = abs(S2_f[t])
    }
    
  }
  
  ####### calculating S1, S2, S3 ########
  Z1 = min(S1_f[n0:n1])
  predicted1 = which.min(S1_f[n0:n1]) + n0-1
  
  if(!correct){ # if does not do correction
    Z2 = min(S2_f[n0:n1])
    predicted2 = which.min(S2_f[n0:n1]) + n0-1
  }else{
    Z2 = min(correct_S2_f[n0:n1])
    predicted2 = which.min(correct_S2_f[n0:n1]) + n0-1
  }
  
  Z3 = min(S3_f[n0:n1])
  predicted3 = which.min(S3_f[n0:n1]) + n0-1
  predictedDubey = which.min(D_f[n0:n1]) + n0-1
  
  return(list(Z1 = Z1, location1 = predicted1,  
              Z2 = Z2, location2 = predicted2, 
              Z2 = Z3, location3 = predicted3))
  
}

