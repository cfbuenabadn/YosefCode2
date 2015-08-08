library(RUVSeq)
runRUVg <- function(e, control_vec = NULL,K) {
  if(!is.integer(e)){
    e = 2^(round(log2(e)))
  }
  if(is.null(control_vec)){
    control_vec = 1:dim(e)[1]
  }
  ruv.obj <- RUVg(e,control_vec,k = K)
  return(ruv.obj[[2]])
}