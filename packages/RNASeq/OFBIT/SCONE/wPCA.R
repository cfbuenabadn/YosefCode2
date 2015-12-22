## ===== Weighted PCA and Weighted Moment Functions =====
# Contains functions to compute weighted means, weighted covariance matrices, 
# and weighted PCA
# Author: Michael Cole, Yosef Lab 2015

library(rARPACK,quietly = T) # ARPACK wrapper

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
  
  SD_EPSILON = 1e10 * .Machine$double.eps #~2.2e-6
  
  z = w*(t(t(x) - colSums(w*x)/colSums(w)))
  Z = t(z) %*% z
  W = t(w) %*% w
  out = Z/W
  out[W < SD_EPSILON] = NA
  return(out)
}

## ----- Weighted PCA
# Returns list containing weighted PCA results.
# Columns correspond to data points, rows to features.
# x = data matrix
# w = weight matrix
# nu = number of eigenvectors to compute
# filt = whether columns should be filtered to prevent weightless products
# scale = if true, scaling is performed for PCA on correlation, vs covariance

wPCA = function(x,w,nu = min(dim(x)), filt = F, scale = T){
  
  SD_EPSILON = 1e10 * .Machine$double.eps #~2.2e-6
  
  # No negative weights allowed.
  stopifnot(!any(w < 0))
  
  # Filter out weightless pairs
  if(filt){
    to.pass = rowSums((w %*% t(w)) < SD_EPSILON) == 0
    x = x[to.pass,]
    w = w[to.pass,]
  }
  
  # Weighted Covariance/Correlation Matrix
  cov = wcov(t(x),t(w))
  stopifnot(!any(is.na(cov))) # Stop if covariance matrix is undefined
  if(scale){
    var = diag(cov)
    cor = cov/sqrt(var%*%t(var))
  }else{
    cor = cov
  }
  
  # Eigenvectors of Weighted Correlation Matrix
  eig_obj = eigs(cor,k = nu,which = "LM")
  eigvecs = eig_obj$vectors
  
  # Project centered and scaled(?) data along principal eigenvectors
  # to compute weighted principal components
  if(scale){
    z = (x - rowSums(x*w)/rowSums(w))/sqrt(var)
  }else{
    z = (x - rowSums(x*w)/rowSums(w))
  }
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