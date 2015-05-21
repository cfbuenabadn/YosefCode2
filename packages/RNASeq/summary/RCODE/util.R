cor2  = function(x,y,corr_method = "spearman"){
  c = cor(x,y,method = corr_method)
  Fish = (1/2)*log((1+c)/(1-c))
  z = sqrt((dim(x)[1]-3)/1.06)*Fish
  p = 2*pnorm(-abs(z))
  return(list(r = c, p = p))
}
