## ===== Weighted PCA and Weighted Moment Functions =====
# Contains functions to compute weighted means, weighted covariance matrices, 
# and weighted PCA
# Author: Michael Cole, Yosef Lab 2015

require(rARPACK) # ARPACK wrapper

## ----- Weighted Mean
# Returns vector containing the weighted mean of an input vector
# x = data vector
# w = weight vector
wmean = function(x,w){
  return(sum(x*w)/sum(w))
}

## ----- Weighted Covariance
# Returns matrix containing the weighted covariance matrix of column vectors. 
# Each column is centered prior to product.
# x = data matrix
# w = weight matrix
# Note: Function will return NA if the weight product of any two columns is zero, 
# as covariance in undefined for the pair.
wcov = function(x,w){
  z = w*(t(t(x) - colSums(w*x)/colSums(w)))
  Z = t(z) %*% z
  W = t(w) %*% w
  out = Z/W
  out[W == 0] = NA
  return(Z/W)
}

## ----- Weighted PCA
# Returns list containing weighted PCA results.
# Columns correspond to data points, rows to features.
# x = data matrix
# w = weight matrix
# nu = number of eigenvectors to compute
# filt = whether columns should be filtered to prevent weightless products

wPCA = function(x,w,nu = min(dim(x)), filt = F){
  
  # No negative weights allowed.
  stopifnot(!any(w < 0))
  
  # Filter out weightless pairs
  if(filt){
    to.pass = rowSums((w %*% t(w)) == 0) == 0
    x = x[to.pass,]
    w = w[to.pass,]
  }
  
  # Weighted Covariance/Correlation Matrix
  cov = wcov(t(x),t(w))
  stopifnot(!any(is.na(cov))) # Stop if covariance matrix is undefined
  var = diag(cov)
  cor = cov/sqrt(var%*%t(var))
  
  # Eigenvectors of Weighted Correlation Matrix
  eig_obj = eigs(cor,k = nu,which = "LM")
  eigvecs = eig_obj$vectors
  
  # Project scaled data (weighted moments) along principal eigenvectors
  # to compute weighted principal components
  z = (x - rowSums(x*w)/rowSums(w))/sqrt(var)
  wpc = t(z*w) %*% eigvecs
  
  # Output Results (prcomp convention)
  res = list()
  res$x = wpc
  res$sdev = eig_obj$values
  res$rotation = eig_obj$vectors
  
  # Return filter vector if filter option was chosen
  if(filt){
    res$to.pass = to.pass
  }
  
  return(res)
}