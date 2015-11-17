library(DESeq)
library(EDASeq)
library(edgeR)
library(sva)
library(aroma.light)
library(lme4)


#' Null function
#' @param ei = Numerical matrix. (rows = genes, cols = samples)
#' @return ei
ID_FN = identity

#' Upper-quartile normalization wrapper.
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return Upper-quartile normalized matrix (scaled by sample UQ)
UQ_FN = function(ei){
  eo = betweenLaneNormalization(ei, which="upper", round = FALSE)
  return(eo)
}

#' Upper-quartile normalization applied to positive data.
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return Upper-quartile (positive) normalized matrix (scaled by sample UQ)
UQ_FN_POS = function(ei){
  zei = ei
  zei[ei == 0] = NA
  q = apply(zei, 2, quantile, 0.75, na.rm = T)
  zeo = t(t(zei)/q)*mean(q)
  eo = zeo
  eo[ei == 0] = 0
  return(eo)
}

#' Full-Quantile normalization wrapper.
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return FQ normalized matrix.
FQ_FN = function(ei){
  eo = normalizeQuantileRank.matrix(ei)
  return(eo)
}

#' Full-Quantile normalization applied to positive data.
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return FQ (positive) normalized matrix.
FQ_FN_POS = function(ei){
  
  # Vector of integers used for computation
  base_rank = cumsum(rep(1,dim(ei)[1]))
  
  # Quantile Index Matrix: Values between 0 and 1 corresponding to quantile
  quant_mat = NULL
  
  # Re-ordered Data Matrix
  x_mat = NULL
  print("Sorting Matrix...")
  # For each sample:
  for (i in 1:dim(ei)[2]){
    # Sort data and replace zeroes with NA
    x_mat = cbind(x_mat,rev(sort(ei[,i])))
    x_mat[x_mat == 0 ] = NA
    # Compute corresponding quantile indices for that sample
    quant = base_rank/sum(ei[,i]>0)
    quant[quant > 1] = NA
    quant_mat = cbind(quant_mat,quant)
  }
  # Vector form of quantile index matrix
  quant_out = as.numeric(as.vector(quant_mat))
  print("Complete.")
  print("Spline Interpolation...")
  # Interpolation Matrix (Values of all quantiles)
  inter_mat = rep(0,length(quant_out))
  ob_counts = rep(0,length(quant_out)) # Number of observations for averaging
  # For each sample
  for (i in 1:dim(ei)[2]){
    # Produce spline interpolation for that sample, value ~ quantile index
    x1 = na.omit(quant_mat[,i])
    y1 = na.omit(x_mat[,i])
    # Evaluated at all possible quantile indices
    inter = approx(x1,y1,xout = quant_out, rule = 2)$y
    ob_counts = ob_counts + !is.na(inter)
    inter[is.na(inter)] = 0
    inter_mat = inter_mat + inter
  }
  print("Complete.")
  
  # Average over the interpolated values from all samples
  inter_mean = inter_mat/ob_counts
  
  ## Substituting Mean Interpolated Values for Expression Values and Return
  eo = matrix(inter_mean,ncol = dim(ei)[2])
  eo[is.na(eo)] = 0
  for (i in 1:dim(ei)[2]){
    eo[,i] = rev(eo[,i])[order(order(ei[,i]))]
  }
  print("Finished!")
  return(eo)
}

#' DESeq normalization wrapper.
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return DESeq normalized matrix (scaled by sample )
DESEQ_FN = function(ei){
  size_fac = estimateSizeFactorsForMatrix(ei)
  eo = t(t(ei)/size_fac)
  return(eo)
}

#' DESeq normalization applied to positive data.
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return DESeq (positive) normalized matrix (scaled by sample )
DESEQ_FN_POS = function(ei){
  if(any(ei < 0)){
    stop("Negative values in input.")
  }
  if(!is.null(dim(ei))){
    y = ei
    y[y == 0] = NA # Matrix with zeroes replaced w/ NA
    geom_mean = exp(apply(log(y),1,sum,na.rm = T)/rowSums(!is.na(y))) # Compute Geometric Mean of Expression for Each Gene (Use positive data only)
  }else{
    stop("Null imput matrix dimension.")
  }
  if(!any(geom_mean > 0)){
    stop("Geometric mean non-positive for all genes.")
  } 
  ratios = ei / geom_mean # Divide each Expression Value by Geometric Mean of Corresponding Gene
  if(any(is.infinite(geom_mean))){
    stop("Infinite mean! This should never happen :-<")
  }
  ratios = ratios[geom_mean > 0,] # Ignore Genes with Zero Mean
  y = ratios
  y[y == 0] = NA
  size = apply(y,2,median,na.rm = T) # Size factor taken as median of ratios (positive data only)
  if(any(size == 0)){
    stop("Zero library size. This should never happen :-(")
  }
  eo = t(t(ei)/size)*mean(size)
  return(eo)
}

#' TMM normalization wrapper.
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @return TMM normalized matrix (scaled by sample )
TMM_FN = function(ei){
  size_fac = calcNormFactors(ei,method = "TMM")
  eo = t(t(ei)/size_fac)
  return(eo)
}

#' Output residuals + intercept of fit to quality score matrix
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @param default_args_list = list of default arguments
#' @param args_list = list of user-provided arguments
#' @return Matrix adjusted for K quality scores

QUALREG_FN = function(ei, default_args_list, args_list){
  EPSILON = default_args_list$EPSILON
  scores = prcomp(default_args_list$q,center = TRUE,scale = TRUE)$x[,1:default_args_list$K]
  re = log(ei + EPSILON)
  for (g in 1:dim(ei)[1]){
    lin.mod = lm(re[g,] ~ scores)
    re[g,] = lin.mod$coefficients[1] + lin.mod$residuals
  }
  re = exp(re)
  re[e == 0] = 0 
  return(re)
}

#' Output residuals + intercept of fit to control-gene derived factors
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @param default_args_list = list of default arguments
#' @param args_list = list of user-provided arguments
#' @return Matrix adjusted for K unwanted factors

RUVG_FN = function(ei, default_args_list, args_list){
  EPSILON = default_args_list$EPSILON
  scores = prcomp(t(log(ei[default_args_list$control_genes,]+default_args_list$EPSILON)),center = TRUE,scale. = FALSE)$x[,1:default_args_list$K]   
  re = log(ei + default_args_list$EPSILON)
  for (g in 1:dim(ei)[1]){
    lin.mod = lm(re[g,] ~ scores)
    re[g,] = lin.mod$coefficients[1] + lin.mod$residuals
  }
  re = exp(re)
  re[e == 0] = 0 
  return(re)
}

#' ComBat Wrapper
#' @param ei = Numerical matrix. (rows = genes, cols = samples). Unique row.names are required.
#' @param default_args_list = list of default arguments
#' @param args_list = list of user-provided arguments
#' @return Matrix adjusted for batch

COMBAT_FN = function(ei, default_args_list, args_list){
  EPSILON = default_args_list$EPSILON
  SD_EPSILON = default_args_list$SD_EPSILON
  ei = log(ei + default_args_list$EPSILON)
  pheno = as.data.frame(cbind(default_args_list$tech_batch,default_args_list$bio_cond))
  if(is.null(default_args_list$bio_cond) || any(is.na(default_args_list$bio_cond))){
    colnames(pheno) = c("batch")
    is_var = apply(ei,1,sd) > SD_EPSILON
    mod = model.matrix(~ 1, data = pheno)
  }else{
    colnames(pheno) = c("batch","phenotype")
    is_var = T
    for (p in unique(default_args_list$bio_cond)){
      is_var = is_var & (apply(ei[,default_args_list$bio_cond == p],1,sd) > SD_EPSILON)
    }
    mod = model.matrix(~ as.factor(phenotype), data = pheno)
  }
  eo = ei
  eo[is_var,] = ComBat(dat=ei[is_var,], batch=default_args_list$tech_batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
  eo = exp(eo)
  eo[ei == 0] = 0
  return(eo)
}

NESTED_BATCH_FIXED_FN = function(ei,bio,batch,uv = NULL,w = NULL){
  if(!is.null(bio)){
    bio = factor(bio, levels = sort(levels(bio)))
    batch = factor(batch, levels = unique(batch[order(bio)]))
  }
  
  if(!is.null(bio)){
    n_vec = NULL
    for(lev in levels(bio)){
      n_vec = c(n_vec,nlevels(droplevels(batch[bio==lev])))
    }
  }else{
    n_vec = nlevels(droplevels(batch))
  }
  
  mat = matrix(0,nrow = sum(n_vec),ncol = sum(n_vec - 1))
  xi = 1
  yi = 1
  for(i in 1:length(n_vec)){
    if(n_vec[i] > 1){
      cs = contr.sum(n_vec[i])
      dd = dim(cs)
      mat[xi:(xi + dd[1] - 1),yi:(yi + dd[2] - 1)] = cs
      xi = xi + dd[1]
      yi = yi + dd[2]
    }else{
      xi = xi + 1
    }
  }
  
  leo = log(ei+1)
  if(is.null(uv)){
    for(i in 1:dim(ei)[1]){
      if(!is.null(bio)){
        lm_object = lm(leo[i,] ~ bio + batch,
                       contrasts=list(bio=contr.sum, batch=mat),
                       weights = w[i,])
      }else{
        lm_object = lm(leo[i,] ~ batch,
                       contrasts=list( batch=mat),
                       weights = w[i,])
      }
      bind = (length(lm_object$coefficients) - (sum(n_vec - 1) - 1)):(length(lm_object$coefficients))
      leo[i,] = leo[i,] - (lm_object$contrasts$batch %*% as.matrix(lm_object$coefficients[bind]))[batch]
    }
  }else{
    for(i in 1:dim(ei)[1]){
      if(!is.null(bio)){
        lm_object = lm(leo[i,] ~ bio + batch + uv,
                       contrasts=list(bio=contr.sum, batch=mat),
                       weights = w[i,])
      }else{
        lm_object = lm(leo[i,] ~ batch + uv,
                       contrasts=list(batch=mat),
                       weights = w[i,])
      }
      uvind = (length(lm_object$coefficients) - dim(uv)[2] + 1):(length(lm_object$coefficients))
      bind = (length(lm_object$coefficients) - dim(uv)[2] - (sum(n_vec - 1) - 1)):(length(lm_object$coefficients) - dim(uv)[2])
      leo[i,] = leo[i,] - (lm_object$contrasts$batch %*% as.matrix(lm_object$coefficients[bind]))[batch] - uv %*% lm_object$coefficients[uvind]
    }
  }
  return(leo)
}

NESTED_BATCH_RANDOM_FN = function(ei,bio,batch,uv = NULL,w = NULL){
  design <- model.matrix(~batch - 1)
  leo = log(ei+1)
  if(is.null(uv)){
    for(i in 1:dim(ei)[1]){
      lm_object = lmer(leo[i,] ~ 1 + bio + (1 | batch),weights = w[i,])
      betaj <- ranef(lm_object)[[1]][,1]
      leo[i,] = leo[i,] - as.vector(design %*% betaj)
    }
  }else{
    for(i in 1:dim(ei)[1]){
      lm_object = lmer(leo[i,] ~ 1 + bio + uv + (1 | batch),
                       weights = w[i,])
      betaj <- ranef(lm_object)[[1]][,1]
      uvind = (length(coef(lm_object)[[1]][1,]) - dim(uv)[2] + 1):(length(coef(lm_object)[[1]][1,]))
      leo[i,] = leo[i,] - as.vector(design %*% betaj) - as.vector(uv %*% t(coef(lm_object)[[1]][1,][uvind]))
    }
  }
  return(leo)
}

FACTORIAL_BATCH_FN = function(ei,bio,batch,uv = NULL,w = NULL){
  contrasts <- contr.treatment(length(levels(batch)))
  rownames(contrasts) = levels(batch)
  leo = log(ei+1)
  
  if(is.null(bio)){
    if(is.null(uv)){
      
      # Just Batch
      for(i in 1:dim(ei)[1]){
        lm_object = lm(leo[i,] ~ batch,
                       weights = w[i,])
        bind = (length(lm_object$coefficients) + (2 - (length(levels(batch))))):(length(lm_object$coefficients))
        leo[i,] = leo[i,] - (contrasts %*% as.matrix(lm_object$coefficients[bind]))[batch]
      }
    }else{
      
      # Batch + UV
      for(i in 1:dim(ei)[1]){
        lm_object = lm(leo[i,] ~ batch + uv,
                       weights = w[i,])
        uvind = (length(lm_object$coefficients) - dim(uv)[2] + 1):(length(lm_object$coefficients))
        bind = (length(lm_object$coefficients) - dim(uv)[2] + (2 - (length(levels(batch))))):(length(lm_object$coefficients) - dim(uv)[2])
        leo[i,] = leo[i,] - (contrasts %*% as.matrix(lm_object$coefficients[bind]))[batch] - uv %*% lm_object$coefficients[uvind]
      }
    }
  }else{
    if(is.null(uv)){
      
      # Batch + BIO
      for(i in 1:dim(ei)[1]){
        lm_object = lm(leo[i,] ~ bio + batch,
                       weights = w[i,])
        bind = (length(lm_object$coefficients) + (2 - (length(levels(batch))))):(length(lm_object$coefficients))
        leo[i,] = leo[i,] - (contrasts %*% as.matrix(lm_object$coefficients[bind]))[batch]
      }
    }else{
      
      # Batch + UV + BIO
      for(i in 1:dim(ei)[1]){
        lm_object = lm(leo[i,] ~ bio + batch + uv,
                       weights = w[i,])
        uvind = (length(lm_object$coefficients) - dim(uv)[2] + 1):(length(lm_object$coefficients))
        bind = (length(lm_object$coefficients) - dim(uv)[2] + (2 - (length(levels(batch))))):(length(lm_object$coefficients) - dim(uv)[2])
        leo[i,] = leo[i,] - (contrasts %*% as.matrix(lm_object$coefficients[bind]))[batch] - uv %*% lm_object$coefficients[uvind]
      }
    }
  }
  return(leo)
}

BATCHFREE_FN = function(ei,bio,uv = NULL,w = NULL){
  leo = log(ei+1)
  
  if(is.null(bio)){
    # Just UV
    for(i in 1:dim(ei)[1]){
      lm_object = lm(leo[i,] ~ uv,
                     weights = w[i,])
      uvind = (length(lm_object$coefficients) - dim(uv)[2] + 1):(length(lm_object$coefficients))
      leo[i,] = leo[i,] - uv %*% lm_object$coefficients[uvind]
    }
  }else{
    # Bio + UV
    for(i in 1:dim(ei)[1]){
      lm_object = lm(leo[i,] ~ bio + uv,
                     weights = w[i,])
      uvind = (length(lm_object$coefficients) - dim(uv)[2] + 1):(length(lm_object$coefficients))
      leo[i,] = leo[i,] - uv %*% lm_object$coefficients[uvind]
    }
  }
  return(leo)
}

ADJUST_FN = function(ei,batch = NULL,
                     bio = NULL,
                     uv = NULL,
                     w = NULL,
                     design = c("factorial","nested"), 
                     nested_model = c("fixed","random")){
  
  if(!is.null(batch)){
    design <- match.arg(design)
    
    if(design == "factorial"){
      leo = FACTORIAL_BATCH_FN(ei,bio,batch,uv,w)
    }else{
      nested_model <- match.arg(nested_model)
      
      if(nested_model == "fixed"){
        leo = NESTED_BATCH_FIXED_FN(ei,bio,batch,uv,w)
      }else{
        leo = NESTED_BATCH_RANDOM_FN(ei,bio,batch,uv,w)
      }
    }
  }else{
    leo = BATCHFREE_FN(ei,bio,uv,w)
  }
  return(leo)
}