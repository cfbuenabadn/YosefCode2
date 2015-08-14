require(rARPACK)

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
wPCA = function(x,w,nu = min(dim(x)),niter = 100){
  
  # Weighted Covariance/Correlation Matrix
  cov = wcov(t(x),t(w))
  var = diag(cov)
  cor = cov/sqrt(var%*%t(var))
  
  # SVD of Weighted Correlation Matrix
  eig_obj = eigs(cor,k = nu,which = "LM")
  eigvecs = eig_obj$vectors
  
  #  Projection (Weighted)
  z = (x - rowSums(x*w)/rowSums(w))/sqrt(var) # Weighted mean-subtracted and scaled to weighted variance
  wpc = t(z*w) %*% eigvecs
  
  # Results
  res = list()
  res$sdev = eig_obj$values
  res$rotation = eig_obj$vectors
  res$x = wpc
  return(res)
  
}