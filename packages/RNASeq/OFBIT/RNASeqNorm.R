source("../summary/RCODE/wPCA.R")
library(sva)
library(preprocessCore)
library(mixtools)
library(diptest)

## ===== SD_EPSILON: Constant for purpose of correlation computation =====
SD_EPSILON = 1e10 * .Machine$double.eps #~2.2e-6

## ----- TFilter: Transcript Filtering Wrapper Function -----
# Perform transcript filtering on expression matrix, producing boolean vector annotating genes that pass filter (T) and fail filter (F)
# Two filtering methods are currently available for count and TPM data ("strong" vs "weak"). "Coverage" has been proposed but not yet implemented.

# e = expression matrix (rows = transcripts, cols = samples)
# cov = coverage data (not used)
# type = data type
# method = filtering method

TFilter = function(e,cov = NULL,
                   type = c("count","TPM"), 
                   method = c("strong","weak","coverage")){
  if(any(is.na(e))){
    stop("Expression matrix contains missing values!")
  }
  
  # Select appropriate thresholds based on options.
  
  STHRESH = NULL # Strong: Inclusive transcript failure threshold (Counts or TPM)
  SPROPFAIL = NULL # Strong: Minimum fraction of samples with failed transcript to filter-out that transcript
  WTHRESH = NULL # Weak: Inclusive transcript failure threshold (Counts or TPM)
  WPROPFAIL = NULL # Weak: Minimum fraction of samples with failed transcript to filter-out that transcript
  type <- match.arg(type)
  if(type == "count"){
    STHRESH = 20
    SPROPFAIL = 0.85
    WTHRESH = 10
    WPROPFAIL = 0.9
  }else{
    STHRESH = 10
    SPROPFAIL = 0.20
    WTHRESH = 1
    WPROPFAIL = 0.10
  }
  
  # Perform filtering: Gene fails if expression value falls below THRESH in more than (inclusive) PROPFAIL of samples
  tf.vec = NULL # Transcript Filter Vector
  method <- match.arg(method)
  if(method == "strong"){
    tf.vec = rowMeans(e <= STHRESH) < SPROPFAIL
  }else if(method == "weak"){
    tf.vec = rowMeans(e <= WTHRESH) < WPROPFAIL
  }else {
    stop("To Be Implemented: Coverage-Based Gene Filtering")
  }
  
  return(tf.vec)
}

## ----- FastComBat: ComBat Wrapper Function -----
# Run empirical Bayes method for batch-adjustment using simple arguments.

# e = expression matrix (rows = transcripts, cols = samples)
# batch = categorical batch variable for adjustment
# biobatch = categorical biological covariate for use in adjustement model
# par.prior = use parametric prior? ComBat argument
# prior.plots = plot priors? ComBat argument
# to.log = should data be logged prior to adjustment? Recommended for expression data.

# WARNING! ComBat will adjust zero data to non-zero values

FastComBat = function(e,batch,biobatch = NULL,par.prior=T,prior.plots=F,to.log = T){
  # Log-Transform Data
  if(to.log){
    e = log(e + 1)
  }
  
  # Pheno Data
  pheno = as.data.frame(cbind(batch,biobatch))
  
  if(is.null(biobatch) || any(is.na(biobatch))){
    # If no biobatch or some samples are not catgorized, use a constant model
    colnames(pheno) = c("batch")
    is_var = apply(e,1,sd) > SD_EPSILON
    mod = model.matrix(~ 1, data = pheno)
  }else{
    # Otherwise, expression is a function of categorical biology
    colnames(pheno) = c("batch","phenotype")
    # Expression must be variable in each biological condition, otherwise inference is problematic
    is_var = T
    for (p in unique(biobatch)){
      is_var = is_var & (apply(e[,biobatch == p],1,sd) > SD_EPSILON)
    }
    mod = model.matrix(~ as.factor(phenotype), data = pheno)
  }
  combat_edata = e
  combat_edata[is_var,] = ComBat(dat=e[is_var,], batch=batch, mod=mod, par.prior=par.prior, prior.plots=prior.plots)
  # It seems to make sense that we leave !is_var transcripts alone.
  
  # Exp-Transform Data
  if(to.log){
    combat_edata = exp(combat_edata) - 1
    combat_edata[combat_edata < 0] = 0
  }

  return(combat_edata)
}

## ----- Upper Quartile Normalization -----
# Returns upper-quartile normalized matrix
# x = expression matrix (rows = transcripts, cols = samples)
UQ = function(x){
  q = apply(x,2,quantile,na.rm = T)[4,]
  y = t(t(x)/q)*mean(q)
  y[,q == 0] = NA
  return(y)
}

## ----- DESeq Normalization by Size Factors ------
# Returns DESeq-scaled normalized matrix
# x = expression matrix (rows = transcripts, cols = samples)
# singlecell = if True, compute geometric means of positive data only

DESeqNorm = function(x, singlecell = F){
  if(any(x < 0)){
    stop("Negative values in input.")
  }
  
  if(!is.null(dim(x))){
    if(singlecell){
      y = x
      y[y == 0] = NA # Matrix with zeroes replaced w/ NA
      geom_mean = exp(apply(log(y),1,sum,na.rm = T)/rowSums(!is.na(y))) # Compute Geometric Mean of Expression for Each Gene (Use positive data only)
    }else{
      y = x
      y[y == 0] = NA # Matrix with zeroes replaced w/ NA
      geom_mean = exp(apply(log(y),1,sum)/(dim(y)[2])) # Compute Geometric Mean of Expression for Each Gene
      geom_mean[is.na(geom_mean)] = 0
    }
  }else{
    stop("Null imput matrix dimension.")
  }
  
  if(!any(geom_mean > 0)){
    stop("Geometric mean non-positive for all genes.")
  }
  
  ratios = x / geom_mean # Divide each Expression Value by Geometric Mean of Corresponding Gene
  
  if(any(is.infinite(geom_mean))){
    stop("Infinite mean! This should never happen :-<")
  }
  
  ratios = ratios[geom_mean > 0,] # Ignore Genes with Zero Mean
  
  if(singlecell){
    y = ratios
    y[y == 0] = NA
    size = apply(y,2,median,na.rm = T) # Size factor taken as median of ratios (positive data only)
  }else{
    size = apply(ratios,2,median) # Size factor taken as median of ratios
  }
  
  if(any(size == 0)){
    stop("Zero library size. This should never happen :-(")
  }
  
  y = t(t(x)/size)
  return(y)
}

## ----- Full Quantile Normalization with Average-Rank Method -----
# Returns FQ normalized matrix using average rank to handle ties / restoring zeroes
# x = expression matrix (rows = transcripts, cols = samples)
MeanQuant = function(x){
  y = normalize.quantiles(x)
  y[x == 0] = 0
  return(y)
}

## ----- Conditional Quantile Normalization -----
# Full Quantile normalization applied to non-zero data, under assumption of no false-positives
# x = expression matrix (rows = transcripts, cols = samples)

ConQuant = function(x,v = F){
  
  # Vector used for computation
  base_rank = cumsum(rep(1,dim(x)[1]))
  
  ## Computing Sample-Specific Quantiles
  print("0/2: Computing Sample-Specific Quantiles")
  
  # Quantile Index Matrix: Values between 0 and 1 corresponding to quantile
  quant_mat = NULL
  # Re-ordered Data Matrix
  x_mat = NULL
  # For each sample:
  for (i in 1:dim(x)[2]){
    # Sort data and replace zeroes with NA
    x_mat = cbind(x_mat,rev(sort(x[,i])))
    x_mat[x_mat == 0 ] = NA
    # Compute corresponding quantile indices for that sample
    quant = base_rank/sum(x[,i]>0)
    quant[quant > 1] = NA
    quant_mat = cbind(quant_mat,quant)
  }
  # Vector form of quantile index matrix
  quant_out = as.numeric(as.vector(quant_mat))
  
  ## Interpolating Quantiles
  print("1/2: Interpolating Quantiles")
  # Interpolation Matrix
  inter_mat = rep(0,length(quant_out))
  ob_counts = rep(0,length(quant_out))
  # For each sample
  for (i in 1:dim(x)[2]){
    if(v == T){print(i)}
    # Produce spline interpolation for that sample, value ~ quantile index
    x1 = na.omit(quant_mat[,i])
    y1 = na.omit(x_mat[,i])
    # Evaluated at all possible quantile indices
    inter = approx(x1,y1,xout = quant_out, rule = 2)$y
    ob_counts = ob_counts + !is.na(inter)
    inter[is.na(inter)] = 0
    inter_mat = inter_mat + inter
  }
  # Average over the interpolated values from all samples
  inter_mean = inter_mat/ob_counts
  
  ## Substituting Mean Interpolated Values for Expression Values and Return
  print("2/2: Substituting Expression Values")
  inter_mat = matrix(inter_mean,ncol = dim(x)[2])
  inter_mat[is.na(inter_mat)] = 0
  for (i in 1:dim(x)[2]){
    inter_mat[,i] = rev(inter_mat[,i])[order(order(x[,i]))]
  }
  return(inter_mat)
}

## ----- Extend Adjustment Ratios from Subset of Genes to All Genes -----
# Extends normalization scheme from a normalized subset of genes to all other genes, 
# keeping zeroes fixed. An un-normalized gene will be assigned the mean ratio of the 
# next lowest rank and next highest rank normalized gene. If only one exists, that ratio
# is assigned.
#
# raw = vector of un-normalized data
# norm = vector of normalized values, with NA for all missing normalized data
# simple = assume the normalized data to be scaled by a constant factor

ExtendRat = function(raw, norm, simple = FALSE){
  
  if (any(is.na(raw))){
    stop("No missing data allowed in raw data input.")
  }
  
  # Simple case
  if(simple){
    rat = na.omit(norm[raw > 0]/raw[raw > 0])[1]
    return(raw*rat)
  }
  
  o = order(raw) # Order of raw data 
  r = norm[o]/raw[o] # Ordered Normalization Ratio Vector
  is_zero = (raw[o] == 0) # Data that is zero in raw
  r = r[!is_zero] # Remove Zeroes (these sit on the left in the ordered vector)
  len = length(r) # Number of ratios for non-zero data values
  
  # Leading NA propagates as -1 (No information will be gained from -1 entries)
  if(is.na(r[1])){
    r[1] = -1
  }
  
  # Left Ratio - For all data, define the left adjustment ratio for missing data as the nearest adjustment ratio on the left
  left.r = r[2:len] # Only consider ratios that could be NA after modifying leading value - remember that r[1] would have been converted to -1 if it was NA!
  na.index = which(is.na(left.r)) # Track whether NA/undefined ratios still persist
  while(length(na.index) > 0){
    left.r[na.index] = c(r[1],left.r)[na.index] # Replace undefined ratio with ratio on the left
    na.index = which(is.na(left.r)) # Update undefined left ratios
  }
  # Produce final left ratio
  left.r = c(r[1],left.r)  # Restore leading value to the left ratio vector
  left.r[left.r == -1] = NA # If it is -1, restore it's value to NA
  
  # Restore Leading NA
  if(r[1] == -1){
    r[1] = NA
  }
  
  ## Same thing, now on the right!
  # Tail NA propagates as -1 (No Information Gained From this Entry)
  if(is.na(r[len])){
    r[len] = -1
  }
  # Right Ratio
  right.r = r[1:(len-1)] # Ratios that could be NA - remember that r[len] would have been converted to -1 if it was NA
  na.index = which(is.na(right.r)) # Track NA/undefined right ratios (missing data)
  while(length(na.index) > 0){
    right.r[na.index] = c(right.r,r[len])[na.index+1] # Replace undefined ratio with ratio on the right
    na.index = which(is.na(right.r)) # Update undefined right ratios
  }
  # Produce final right ratio
  right.r = c(right.r,r[len])
  right.r[right.r == -1] = NA
  # Compute final adjustment ratio as mean of left and right, ignoring NA, and add zeroes back in!
  new.r = c(rep(0,sum(is_zero)),rowMeans(cbind(left.r,right.r),na.rm = T)) 
  # Apply ratio to raw data to produce final normalized data - note that zeroes MUST have been restored.
  return((new.r*raw[o])[order(o)])
} 

## ----- Expression Scaling Wrapper -----
# Perform expression scaling using one of three different scaling (global or quantile-specific) methods, 
# including or excluding zeroes. Normalization scheme is based on tf.vec positive transcripts, extended
# to all others according to similarity in rank.
#
# e = expression matrix (rows = transcripts, cols = samples)
# tf.vec = transcript filter vector
# method = scaling method
# zero.method = zero-handling method

EScale = function(e, tf.vec = T, method = c("FQ","UQ","DESeq"), zero.method = c("all","positive")){
  method <- match.arg(method)
  if(method == "FQ"){
    zero.method <- match.arg(zero.method)
    if(zero.method == "all"){
      no = MeanQuant(e[tf.vec,])
    }else{
      no = ConQuant(e[tf.vec,])
    } 
  }else if(method == "UQ"){
    zero.method <- match.arg(zero.method)
    if(zero.method == "all"){
      no = UQ(e[tf.vec,])
    }else{
      # Ignore zeroes in computing upper quantile
      de = e[tf.vec,]
      de[de == 0] = NA
      no = UQ(de)
      no[de == 0] = 0
    } 
  }else {
    # DESeq Size Factor
    no = DESeqNorm(e[tf.vec,])
  }
  
  # Extend Ratios from Reference Genes
  if(any(!tf.vec)){
    print("Extending Ratios from Reference Genes")
    fno = matrix(NA,nrow = dim(e)[1],ncol = dim(e)[2])
    fno[tf.vec,] = no
    for (s in 1:dim(e)[2]){
      fno[,s] = ExtendRat(raw = e[,s],norm = fno[,s],simple = (method != "FQ"))
    }
    no = fno
  }
  return(no)
}

## ----- Standardize Quality Matrix -----
# Standardize quality metric
# q = quality metric matrix (columns = named features, rows = samples)
# ... = lists for specific transformations (see below)

PPQual = function(q, to.log = c("NREADS", "NALIGNED"), 
                  to.abs.log = c("MEDIAN_5PRIME_TO_3PRIME_BIAS","MEDIAN_5PRIME_BIAS","MEDIAN_3PRIME_BIAS"),
                  to.logit.one = c("PCT_RIBOSOMAL_BASES","PCT_CODING_BASES","PCT_UTR_BASES",
                                   "PCT_INTRONIC_BASES","PCT_INTERGENIC_BASES","PCT_MRNA_BASES"),
                  to.logit.hund = c("RALIGN")){
  
  if(any(is.null(colnames(q)))){
    stop("No quality parameter names.")
  }
  
  quality.features = q
  
  # Convert non-numeric data fields
  for (i in 1:dim(quality.features)[2]){
    if(!is.numeric(quality.features[,i])){
      quality.features[,i] = as.numeric(as.character(quality.features[,i]))
    }
  }
  
  ## ----- Special Transformations
  # 0 to Infinity -> log
  quality.features[,to.log] = log(quality.features[,to.log]+.01)
  # 0 to Infinity, Best at 1
  quality.features[,to.abs.log] = exp(abs(log(quality.features[,to.abs.log]+.01)))-.01
  # 0 to 1
  quality.features[,to.logit.one] = log((quality.features[,to.logit.one]+.01)/(1-quality.features[,to.logit.one]+.01))
  # 0 to 100
  quality.features[,to.logit.hund] = log((quality.features[,to.logit.hund]+1)/(100-quality.features[,to.logit.hund] + 1))
  
  ## ----- Remove NA (missing data), Constants, and scale
  quality.features = t(na.omit(t(quality.features)))
  quality.features = quality.features[,apply(quality.features,2,sd) > SD_EPSILON]
  quality.features = scale(quality.features,center = T,scale = T)
  
  return(quality.features)
}

## ----- Select Processed Quality Features showing Significant Association w/ Expression PCA-----
# Return feature names for implicated quality metrics

# e = expression matrix (rows = transcripts, cols = samples)
# ppq = preprocessed quality metric matrix (continous!)
# tf.vec = transcript filter vector - only transcripts passing filter included in expression PCA
# MAX_EXP_PCS = max number of expression PCs tested for correlations
# qval_thresh = FDR q-value threshold used for assigning significant association between expression PCs and quality metrics
# method = correlation computation method, "spearman" or "pearson"

QSel = function(e, ppq, tf.vec = T, MAX_EXP_PCS = 5, qval_thresh = .01, method = c("spearman","pearson")){
  method <- match.arg(method)
  
  tf.vec = tf.vec & (apply(e,1,sd) > SD_EPSILON) # Only consider variable genes for PCA
  epc = prcomp(t(log(e[tf.vec,] + 1)), center = TRUE, scale = TRUE) # ePCA
  epc_q_cor = cor(epc$x[,1:MAX_EXP_PCS],ppq, method = method) # ePCA vs ppq
  
  # Fisher Transformation and P-value assignment
  fish = (1/2)*log((1+epc_q_cor)/(1-epc_q_cor))
  if(method == "spearman"){
    z = sqrt((dim(epc$x)[1]-3)/1.06)*fish
  }else{
    z = sqrt(dim(epc$x)[1]-3)*fish
  }
  p = 2*pnorm(-abs(z))
  
  # Keep only features showing significant association
  to.keep = p.adjust(p,method = "BH") < qval_thresh
  to.keep.vec = colSums(matrix(to.keep,nrow = MAX_EXP_PCS)) > 0
  return(colnames(ppq)[to.keep.vec])
}

## ----- Generate Quality PCA Scores -----
# Produce prcomp object containg qPCA Scores 

# e = expression matrix (rows = transcripts, cols = samples)
# q = quality metric matrix (columns = named features, rows = samples)
# sig.test = whether or not significance tests should be used to refine list of quality features
# tf.vec = transcript filter vector - only transcripts passing filter included in quality feature selection

QPCScores = function(e,q,sig.test = T,tf.vec = T){
  ppq = PPQual(q)
  
  feat_names = NULL
  if(sig.test){
    feat_names = QSel(e = e,ppq = ppq)
  }else{
    feat_names = colnames(ppq)
  }
  qpc = prcomp(ppq[,feat_names],center = T, scale = T)
  return(qpc)
}

## ----- Generate Quality Score Bin -----
# Generate discretized quality covariates from a matrix of quality scores, using one of two methods,
# quantiles or normal mixture

# score = numeric matrix of quality scores (or others...), columns corresponding to scores
# method = bin method
# mix_k = number of mixture components
# quantile_bins = number of quantile-based bins (equally binned... roughly)

QBins = function(score,method = c("quantiles","mix_norm"), 
                 mix_k = 2, 
                 quantile_bins = 2,
                 dip_thresh = .01){
  
  method <- match.arg(method)
  if(method == "mix_norm" ){
    if(dip.test(score)$p.value < dip_thresh){
      # Membership to Mixture Components 
      nmem = normalmixEM(score,k = mix_k,)
      mp = apply(nmem$posterior,1,max)
      memb = NULL
      for (i in 1:length(mp)){
        memb[i] = which(nmem$posterior[i,] == mp[i]) 
      }
      return(memb)
    }else{
      return(rep(1,length(score)))
    }
  }else{
    # Membership to Partitions
    memb = cut(rank(score),quantile_bins,labels = F)
    return(memb)
  }
}

## ===== Adjustment Methods =====

## ----- Local-Scaling Normalization (Genes) -----
# FOR EACH GENE: Scale the values of each batch by the median, mean or UQ of that batch, 
# relative to all samples. All summary values are produced from un-logged data!

# e = expression matrix (rows = transcripts, cols = samples)
# batch = categorical batch variable
# method = summary value for scaling
# EPSILON = pseudocount for internal log-transformations
# zero.method = zero-handling method: Ignore zeroes? If so, pseudp-count is set to zero.

LocScale = function(e,batch, method = c("median","mean","UQ"), EPSILON = 1, zero.method = c("all","positive")){
  zero.method <- match.arg(zero.method)
  
  if(zero.method == "positive"){
    e[e == 0] = NA
    EPSILON = 0
  }
  
  method <- match.arg(method)
  batches = unique(batch)
  if(method == "median"){
    global.scales = apply(e,1,median, na.rm = T)
  }else if(method == "mean"){
    global.scales = apply(e,1,mean, na.rm = T)
  }else{
    global.scales = apply(e,1,quantile, na.rm = T)[4,]
  }
  ce = e
  for (b in batches){
    bin.samples = (batch == b)
    if(method == "median"){
      if(sum(bin.samples) > 1){
        local.scales = apply(e[,bin.samples],1,median, na.rm = T)
      }else{
        local.scales = median(e[,bin.samples], na.rm = T)
      }
      
    }else if(method == "mean"){
      if(sum(bin.samples) > 1){
        local.scales = apply(e[,bin.samples],1,mean, na.rm = T)
      }else{
        local.scales = mean(e[,bin.samples], na.rm = T)
      }
      
    }else{
      if(sum(bin.samples) > 1){
        local.scales = apply(e[,bin.samples],1,quantile, na.rm = T)[4,]
      }else{
        local.scales = quantile(e[,bin.samples], na.rm = T)[4,]
      }
    }
    new_values = log(e[,bin.samples]+EPSILON) + log(global.scales + EPSILON) - log(local.scales + EPSILON)
    ce[,bin.samples] = exp(new_values) - EPSILON 
    ce[,bin.samples][ce[,bin.samples] < 0] = 0
  } 
  
  if(zero.method == "positive"){
    ce[e == 0] = 0
  }
  
  return(ce)
}

## ----- Residuals of Log-Linear Fit to Genes -----
# Output residuals + intercept of fit to score matrix

# e = expression matrix (rows = transcripts, cols = samples)
# score = numeric matrix of quality scores (or others...), columns corresponding to scores
# EPSILON = pseudocount for internal log-transformations
# zero.method = zero-handling method: Ignore zeroes? If so, pseudp-count is set to zero.

ResLoc = function(e,scores, EPSILON = 1, zero.method = c("all","positive")){
  zero.method <- match.arg(zero.method)
  
  if(zero.method == "positive"){
    w = (e > 0) + 0
    EPSILON = 0
  }else{
    w = matrix(1,nrow = dim(e)[1],ncol = dim(e)[2])
  }

  re = log(e + EPSILON)
  for (g in 1:dim(e)[1]){
    if(sum(w[g,]) > 0){
      lin.mod = lm(re[g,] ~ scores,weights = w[g,])
      re[g,] = lin.mod$coefficients[1] + lin.mod$residuals
    }
  }
  

  re = exp(re) - EPSILON
  re[re < 0] = 0
  
  if(zero.method == "positive"){
    re[e == 0] = 0
  }
  
  return(re)
}

# ----- YL_RUVg -----
# Output residuals + intercept of fit to W matrix produced via RUVg

# control_vec = boolean vector (T = control genes unaffected by biological factors)
# K = number of hidden unwanted factors
# score = numeric matrix of quality scores (or others...), columns corresponding to scores
# zero.method = zero-handling method: Ignore zeroes? If so, pseudp-count is set to zero.

YL_RUVg <- function(e, control_vec = NULL,K, zero.method = c("all","positive")) {
  
  # If no control genes are provided... use all genes!
  if(is.null(control_vec)){
    control_vec = 1:dim(e)[1]
  }
  
  # Derive hidden factors from control genes
  cmat = e[control_vec,]
  zero.method <- match.arg(zero.method)
  if(zero.method == "positive"){
    # Hidden factors in non-zero space
    W = wPCA(log(cmat+1),0 + (cmat > 0),filt = T,nu = K,scale. = F)$x
  }else{
    # Hidden factors in full space
    W = prcomp(t(log(cmat+1)),center = T,scale. = F)$x[,1:K]
  }
  
  return(ResLoc(e = e,scores = W,zero.method = zero.method ))
}