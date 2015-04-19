# Weighted mean
wmean = function(x,w){
  return(sum(x*w)/sum(w))
}

# Weighted Covariance
wcov = function(x,w){
  z = w*(t(t(x) - colSums(w*x)/colSums(w)))
  Z = t(z) %*% z
  W = t(w) %*% w
  return(Z/W)
}

# Weighted Covariance
wcor = function(x,w){
  cov = wcov(t(x),t(w))
  var = diag(cov)
  cor = cov/sqrt(var%*%t(var))
  return(cor)
}

# Weighted PCA
# Rows = features
# Columns = experiments
wPCA = function(x,w){
  
  # Weighted Covariance/Correlation Matrix
  cov = wcov(t(x),t(w))
  var = diag(cov)
  cor = cov/sqrt(var%*%t(var))
  
  # SVD of Weighted Correlation Matrix
  svd_obj = svd(cor)
  eigvecs = svd_obj$u # SVD Takes a while!
  
  #  Projection (Weighted)
  z = (x - rowSums(x*w)/rowSums(w))/sqrt(var) # Weighted mean-subtracted and scaled to weighted variance
  wpc = t(z*w) %*% eigvecs
  
  # Results
  res = list()
  res$sdev = svd_obj$d
  res$rotation = svd_obj$u
  res$x = wpc
  return(res)
  
}