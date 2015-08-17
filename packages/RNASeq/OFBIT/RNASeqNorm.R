## 2) ===== Transcript Filtering =====
TFilter = function(e,cov = NULL,
                   type = c("count","TPM"), 
                   method = c("strong","weak","coverage")){
  if(any(is.na(e))){
    stop("Expression matrix contains missing values!")
  }
  STHRESH = NULL # Strong: Inclusive transcript failure threshold (CPM or TPM)
  SPROPFAIL = NULL # Strong: Minimum fraction of samples with failed transcript to filter-out that transcript
  WTHRESH = NULL # Weak: Inclusive transcript failure threshold (CPM or TPM)
  WPROPFAIL = NULL # Weak: Minimum fraction of samples with failed transcript to filter-out that transcript
  type <- match.arg(type)
  if(type == "count"){
    #e = t(t(e)/colSums(e))*(10^6) # Conversion to counts per million counts (CPM)
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
  
  tf.vec = NULL # Transcript Filter Vector: F == failed transcript
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

## 3) ===== ComBat without Covariates =====
FastComBat = function(e,batch,biobatch = NULL,par.prior=T,prior.plots=F,to.log = T,SD_EPSILON = 0){
  if(to.log){
    e = log(e + 1)
  }
  require(sva)
  if(is.null(biobatch)){
    is_var = apply(e,1,sd) > SD_EPSILON
    pheno = as.data.frame(cbind(batch))
    colnames(pheno) = c("batch")
    batch = pheno$batch
    mod = model.matrix(~ 1, data = pheno)
    combat_edata = e
    combat_edata[is_var,] = ComBat(dat=e[is_var,], batch=batch, mod=mod, par.prior=par.prior, prior.plots=prior.plots)
    combat_edata[!is_var,] = 0
  }else{
    is_var = T
    for (p in unique(biobatch)){
      is_var = is_var & (apply(e[,biobatch == p],1,sd) > SD_EPSILON)
    }
    pheno = as.data.frame(cbind(batch,biobatch))
    colnames(pheno) = c("batch","phenotype")
    batch = pheno$batch
    mod = model.matrix(~ as.factor(phenotype), data = pheno)
    combat_edata = e
    combat_edata[is_var,] = ComBat(dat=e[is_var,], batch=batch, mod=mod, par.prior=par.prior, prior.plots=prior.plots)
    combat_edata[!is_var,] = 0  }
  if(to.log){
    combat_edata = exp(combat_edata) - 1
  }
  return(combat_edata)
}

## ===== Upper Quartile Normalization ====
UQ = function(x, to.log = F ){
  if(to.log){
    x = log(x + 1)
  }
  q = apply(x,2,quantile,na.rm = T)[4,]
  y = t(t(x)/q)*mean(q)
  if(to.log){
    y = exp(y) - 1
  }
  return(y)
}

## ===== DESeq Normalization by Size Factors ====
# Provide Linear-Scale Data!
DESeqNorm = function(x){
  geom_mean = apply(x,1,prod)^(1/dim(x)[2]) # Compute Geometric Mean of Expression for Each Gene
  ratios = x / geom_mean # Divide each Expression Value by Geometric Mean of Corresponding Gene
  ratios = ratios[geom_mean > 0,] # Ignore Genes with Zero Mean
  size = apply(ratios,2,median) # Size factor taken as median of ratios
  y = t(t(x)/size)
  return(y)
}

## ===== Full Quantile Normalization with Average-Rank Method ====
library(preprocessCore)
MeanQuant = function(x, to.log = F ){
  if(to.log){
    x = log(x + 1)
  }
  y = normalize.quantiles(x)
  if(to.log){
    y = exp(y) - 1
  }
  return(y)
}

## ===== Conditional Quantile Normalization ====
# Full Quantile normalization applied to non-zero data, under assumption of no false-positives
ConQuant = function(x,to.log = F){
  if(to.log){
    x = log(x + 1)
  }
  
  base_rank = cumsum(rep(1,dim(x)[1]))
  
  print("0/2: Computing Sample Quantiles")
  quant_mat = NULL
  x_mat = NULL
  for (i in 1:dim(x)[2]){
    x_mat = cbind(x_mat,rev(sort(x[,i])))
    x_mat[x_mat == 0 ] = NA
    
    quant = base_rank/sum(x[,i]>0)
    quant[quant > 1] = NA
    quant_mat = cbind(quant_mat,quant)
  }
  quant_out = as.numeric(as.vector(quant_mat))
  
  print("1/2: Interpolating Quantiles")
  inter_mat = NULL
  for (i in 1:dim(x)[2]){
    x1 = na.omit(quant_mat[,i])
    y1 = na.omit(x_mat[,i])
    inter = approx(x1,y1,xout = quant_out, rule = 2)$y
    inter_mat =cbind(inter_mat,inter)
  }
  inter_mean = rowMeans(inter_mat,na.rm = T)
  
  print("2/2: Substituting Expression Values")
  inter_mat = matrix(inter_mean,ncol = dim(x)[2])
  inter_mat[is.na(inter_mat)] = 0
  for (i in 1:dim(x)[2]){
    inter_mat[,i] = rev(inter_mat[,i])[order(order(x[,i]))]
  }
  
  if(to.log){
    inter_mat = exp(inter_mat) - 1
  }
  return(inter_mat)
}

## ===== Extend Adjustment Ratios from Subset of Genes to All Genes ----
# raw = vector of un-normalized data
# norm = vector of normalized values, with NA for all missing normalized data

ExtendRat = function(raw, norm){
  
  if (any(is.na(raw))){
    stop("No missing data allowed in raw data input.")
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

## 4) ===== Expression Scaling =====
EScale = function(e, fnr = NULL, tf.vec = T, method = c("FQ","UQ","DESeq"), quantile.method = c("mean_rank","positive_only")){
  method <- match.arg(method)
  if(method == "FQ"){
    quantile.method <- match.arg(quantile.method)
    if(quantile.method == "mean_rank"){
      # Average Rank Tie Breaking
      no = MeanQuant(e[tf.vec,],to.log = F)
    }else{
      # Conditional QN - Positive Data Only
      no = ConQuant(e[tf.vec,],to.log = F)
    } 
  }else if(method == "UQ"){
    # Upper Quartile
    no = UQ(e[tf.vec,],to.log = F)
  }else {
    # DESeq Size Factor
    no = DESeqNorm(e[tf.vec,])
  }
  
  # Extend Ratios from Reference Genes
  # TODO: If not full quantile, then simply compute the scaling factor and multiply
  if(any(!tf.vec)){
    print("Extending Ratios from Reference Genes")
    fno = matrix(NA,nrow = dim(e)[1],ncol = dim(e)[2])
    fno[tf.vec,] = no
    for (s in 1:dim(e)[2]){
      print(s)
      fno[,s] = ExtendRat(raw = e[,s],norm = fno[,s])
    }
    no = fno
  }
  return(no)
}

## ===== Standardize Quality Matrix =====
PPQual = function(q, to.log = c("NREADS", "NALIGNED"), 
                  to.abs.log = c("MEDIAN_5PRIME_TO_3PRIME_BIAS","MEDIAN_5PRIME_BIAS","MEDIAN_3PRIME_BIAS"),
                  to.logit.one = c("PCT_RIBOSOMAL_BASES","PCT_CODING_BASES","PCT_UTR_BASES",
                                   "PCT_INTRONIC_BASES","PCT_INTERGENIC_BASES","PCT_MRNA_BASES"),
                  to.logit.hund = c("RALIGN"),
                  SD_EPSILON = 0){
  
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
  
  ## ----- Remove NA, Constants, and scale
  quality.features = t(na.omit(t(quality.features)))
  quality.features = quality.features[,apply(quality.features,2,sd) > SD_EPSILON]
  quality.features = scale(quality.features,center = T,scale = T)
  
  return(quality.features)
}

## ===== Select Processed Quality Features showing Significant Association w/ Expression =====
QSel = function(e, ppq, tf.vec = T, MAX_EXP_PCS = 5, qval_thresh = .01, method = c("spearman","pearson"),SD_EPSILON = 0){
  method <- match.arg(method)
  
  tf.vec = tf.vec & (apply(e,1,sd) > SD_EPSILON) # Only consider variable genes for PCA
  epc = prcomp(t(log(e[tf.vec,] + 1)), center = TRUE, scale = TRUE)
  epc_q_cor = cor(epc$x[,1:MAX_EXP_PCS],ppq, method = method)
  fish = (1/2)*log((1+epc_q_cor)/(1-epc_q_cor))
  if(method == "spearman"){
    z = sqrt((dim(epc$x)[1]-3)/1.06)*fish
  }else{
    z = sqrt(dim(epc$x)[1]-3)*fish
  }
  p = 2*pnorm(-abs(z))
  to.keep = p.adjust(p,method = "BH") < qval_thresh
  to.keep.vec = colSums(matrix(to.keep,nrow = MAX_EXP_PCS)) > 0
  return(colnames(ppq)[to.keep.vec])
}

## 5) ===== Generate Quality PCA Scores =====
QPCScores = function(e,q,sig.test = T,tf.vec = T,center = T, scale = T){
  ppq = PPQual(q)
  
  feat_names = NULL
  if(sig.test){
    feat_names = QSel(e = e,ppq = ppq)
  }else{
    feat_names = colnames(ppq)
  }
  qpc = prcomp(ppq[,feat_names],center = center, scale = scale)
  return(qpc)
}

## ===== Generate Quality Score Bin =====
library(mixtools)
library(diptest)
QBins = function(score,method = c("quantiles","mix_norm","combined"), 
                 mix_k = 2, 
                 quantile_bins = 2, 
                 dip_thresh = .01){
  
  method <- match.arg(method)
  if(method == "mix_norm" || (method == "combined" && (dip.test(score)$p.value < dip_thresh))){
    #print("Membership to Mixture Components")
    nmem = normalmixEM(score,k = mix_k)
    mp = apply(nmem$posterior,1,max)
    memb = NULL
    for (i in 1:length(mp)){
      memb[i] = which(nmem$posterior[i,] == mp[i]) 
    }
    return(memb)
  }else{
    #print("Membership to Partitions")
    memb = cut(rank(score),quantile_bins,labels = F)
    return(memb)
  }
}

## ===== Local-Scaling Normalization (Genes) =====
LocScale = function(e,batch, method = c("median","mean","UQ"), EPSILON = 1){
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
  } 
  return(ce)
}

## ===== Global-Scaling Normalization =====
GlobScale = function(e,batch, method = c("median","mean","UQ"), EPSILON = 1){
  method <- match.arg(method)
  batches = unique(batch)
  if(method == "median"){
    global.scale = median(e)
  }else if(method == "mean"){
    global.scale = mean(e)
  }else{
    global.scale = quantile(e)[4]
  }
  ce = e
  for (b in batches){
    bin.samples = (batch == b)
    if(method == "median"){
      local.scale = median(e[,bin.samples])
    }else if(method == "mean"){
      local.scale = mean(e[,bin.samples])
    }else{
      local.scale = quantile(e[,bin.samples])[4]
    }
    new_values = log(e[,bin.samples]+EPSILON) + log(global.scale + EPSILON) - log(local.scale + EPSILON)
    ce[,bin.samples] = exp(new_values) - EPSILON 
  } 
  return(ce)
}

## ===== Re-Centering Eigen-Genes By Batch =====
AdjEig = function(e,batch, method = c("median","mean"), to.log = T, EPSILON = 1,SD_EPSILON = 0){
  method <- match.arg(method)
  batches = unique(batch)
  
  is_var = apply(e,1,sd) > SD_EPSILON
  oe = e
  if(to.log){
    pc = prcomp(t(log(e[is_var,]+EPSILON)),center = F,scale = F)
  }else{
    pc = prcomp(t(e[is_var,]),center = F,scale = F)
  }
  
  x = t(pc$x)
  if(method == "median"){
    global.centers = apply(x,1,median)
  }else {
    global.centers = apply(x,1,mean)
  }
  
  cx = x
  for (b in batches){
    bin.samples = (batch == b)
    if(method == "median"){
      if(sum(bin.samples) > 1){
        local.centers = apply(x[,bin.samples],1,median)
      }else{
        local.centers = median(x[,bin.samples])
      }
      
    }else{
      if(sum(bin.samples) > 1){
        local.centers = apply(x[,bin.samples],1,mean)
      }else{
        local.centers = mean(x[,bin.samples])
      }
    }
    cx[,bin.samples]= x[,bin.samples] + global.centers - local.centers
  } 
  
  if(to.log){
    ce = exp(pc$rotation %*% cx) - EPSILON
  }else{
    ce = pc$rotation %*% cx
  }
  
  oe[is_var,] = ce
  
  return(oe)
}

## ===== Residuals of Log-Linear Fit to Genes =====
QPCResLoc = function(e,scores, EPSILON = 1, to.log = T){
  if (to.log){
    e = log(e + EPSILON)
  }
  
  re = e
  for (g in 1:dim(e)[1]){
    lin.mod = lm(e[g,] ~ scores)
    re[g,] = lin.mod$residuals
  }
  
  if (to.log){
    re = exp(re) - EPSILON
  }
  return(re)
}

## ===== Residuals of Log-Linear Fit to EigenGenes =====
QPCResEig = function(e,scores, EPSILON = 1, to.log = T,SD_EPSILON = 0){
  if (to.log){
    e = log(e + EPSILON)
  }
  
  oe = e
  is_var = apply(e,1,sd) > SD_EPSILON
  pc = prcomp(t(e[is_var,]),center = F,scale = F)
  oe[!is_var,] = 0
  
  x = t(pc$x)
  rx = x
  for (g in 1:dim(x)[1]){
    lin.mod = lm(x[g,] ~ scores)
    rx[g,] = lin.mod$residuals
  }
  
  oe[is_var,] = pc$rotation %*% rx
  
  if (to.log){
    oe = exp(oe) - EPSILON
  }
  return(oe)
}

## ===== Residuals of Log-Linear Fit to All Genes =====
QPCResGlob = function(e,scores, EPSILON = 1, to.log = T){
  scores = as.matrix(scores)
  if (to.log){
    e = log(e + EPSILON)
  }
  
  re = as.vector(e)
  ascores = NULL
  for (s in 1:dim(scores)[1]){
    mscores = matrix(rep(scores[s,],dim(e)[1]),nrow =dim(e)[1])
    ascores = rbind(ascores,mscores)
  }
  lin.mod = lm(re ~ ascores)
  re = matrix(lin.mod$residuals,nrow = dim(e)[1])
  
  
  if (to.log){
    re = exp(re) - EPSILON
  }
  return(re)
}

# ===== RUVg =====
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

# ===== Metric Normalization =====
library(far)
MetNorm = function(e,scores,to.log = T, K = 50, SD_EPSILON = 0){
  if(to.log){
    e = log(e+1)
  }
  
  # Calculate Distance in Expression
  dev = as.vector(dist(t(e)))
  # Calculate Distances in Score Params
  dqv_tab = NULL
  for(i in 1:dim(scores)[2]){
    dqv = as.vector(dist(scores[,i]))
    dqv_tab = cbind(dqv_tab,dqv)
  }
  
  # Calculate residuals of fit of Distances to Quality Diffs
  lin.mod = lm(dev ~ dqv_tab)
  dev2 = lin.mod$residuals + lin.mod$coefficients[1] 
  # Reconstruct Distance Matrix
  S <- diag(dim(e)[2])
  S = S*0
  S[lower.tri(S, diag=F)] <- dev2
  S = S + t(S)
  d = as.dist(S)
  
  print("Generating Factors...")
  # MDS to generate robust factors (why 3?)
  cmd = cmdscale(d = d,k = dim(S)[1]-1,eig = T)
  cmd = cmd$points[,(cmd$eig[1:dim(cmd$points)[2]] > 10^-10) &  1:dim(cmd$points)[2] <= K]
  ecmd = cmd
  
  print("Generating Eigengenes...")
  
  ## Generate Normalized Expression Matrix
  # Start with eigengenes
  oe = e
  is_var = apply(oe,1,sd) > SD_EPSILON
  pc = prcomp(t(oe[is_var,]),center = F,scale = F)
  x = t(pc$x)
  px = x
  
  print("Fitting Eigengenes to Orthogonal Factors...")
  
  for (g in 1:dim(x)[1]){
    if(g <= dim(cmd)[2]){
      # Replace eigengene with fitted value
      lin.mod = lm(x[g,] ~ ecmd)
      px[g,] = lin.mod$fitted.values 
      # Generate Independent Factors
      if(g < dim(cmd)[2]){
        lin.mod$coefficients[2:(dim(ecmd)[2]+1)]
        v = lin.mod$coefficients[2:(dim(ecmd)[2]+1)]
        v = v/sqrt(sum(v^2))
        rot = diag(dim(ecmd)[2])
        rot[,1] = v
        rot = orthonormalization(rot,norm = T)
        ecmd = (ecmd %*% rot)[,2:(dim(ecmd)[2])]
      }
    }else if(g > dim(cmd)[2]){
      # Ran out of Factors
      px[g,] = 0
    }
  }
  # Rotate back to genes
  oe[is_var,] = pc$rotation %*% px
  # Ground at zero
  oe[oe < 0] = 0
  if(to.log){
    # Restore linear scale
    oe = exp(oe) - 1
  }
  return(oe)
}

## ===== Doubly ZI-Poisson RUVg =====
# e = expression matrix, preferably integer
# control_vec = vector of indexes for control genes
# K = number of hidden factors of unwanted variation
# M = number of hidden drop-out determinants (factors)
# num.iter = number of fit iterations ...
# epsilon = 1

DZIPRUV = function(e,control_vec,K = 1, M = 1, num.iter = 10,epsilon = 1, to.plot = F, impute.only = F){
  data = t(round(e))
  
  print("Estimating hidden factors...")
  ## Finding Hidden Factors of Unwanted Variation
  # Remove Gene-Specific Offsets (Gene Means)
  Z = log(data+epsilon)
  Zs = t(t(Z) - colMeans(Z))
  
  # Find W
  svd_obj = svd(Zs[,control_vec])
  svd_obj$d[1:length(svd_obj$d) > K] = 0
  W = (svd_obj$u %*% diag(svd_obj$d))[,1:K]
  
  ## Finding Hidden Factors of Drop-Out Determination
  # Remove Sample-Specific Offsets
  Z = t(-1 + 2*(data == 0))
  Z = t(t(Z) - colMeans(Z))
  # Find oW
  svd_obj = svd(Z)
  oW = (svd_obj$u)[,1:M]
  
  print("Initializing parameters of poisson component...")
  
  ## Initial Fit For Poisson Component
  O = t(matrix(rep(colMeans(log(data+epsilon)),dim(data)[1]),ncol = dim(data)[1]))
  Z = log(data+epsilon) - O
  ALPHA = NULL
  for(j in 1:dim(data)[2]){
    count = data[,j]
    m1 <- glm(count ~ 0 + W + offset(O[,j]),family = poisson(log))
    coefs = m1$coefficients[1:K]
    coefs[is.na(coefs)] = 0
    ALPHA = cbind(ALPHA,coefs)
  }
  leY = W %*% ALPHA + O
  eY = exp(leY)
  
  print("Initializing parameters of drop-out component...")
  
  ## Initial Fit to Drop Component
  drop_zero_prob = NULL
  for(i in 1:dim(data)[1]){
    success = as.vector((data[i,] == 0))
    failure = as.vector((data[i,] > 0))
    logit_obj = glm(cbind(success,failure) ~ oW,family=quasibinomial(logit))
    drop_zero_prob = rbind(drop_zero_prob,logit_obj$fitted.values)
  }
  
  print("Initializing parameters of off component...")
  
  ## Initial Fit to Off Component
  off_zero_prob = NULL
  BI = NULL
  for(j in 1:dim(data)[2]){
    success = as.vector((data[,j] == 0))
    failure = as.vector((data[,j] > 0))
    logit_obj = glm(cbind(success,failure) ~ 1 ,family=quasibinomial(logit))
    BI[j] = logit_obj$coefficients[1]
    off_zero_prob = cbind(off_zero_prob,logit_obj$fitted.values)
  }
  if(to.plot){
    plot(colMeans(data > 0),1-(1/(1+exp(-BI))))
  }
  
  print(paste("Running",num.iter,"iterations"))
  for (step in 1:num.iter){
    print(paste("Iteration",step,"..."))
    
    # Compute Probability of Zero Assignment
    pois_zero_prob = ppois(0,eY)
    post_drop = drop_zero_prob*(1-off_zero_prob)/((pois_zero_prob*(1-drop_zero_prob) + drop_zero_prob)*(1-off_zero_prob) + off_zero_prob)
    post_drop[data > 0] = 0
    
    post_off = off_zero_prob/((pois_zero_prob*(1-drop_zero_prob) + drop_zero_prob)*(1-off_zero_prob) + off_zero_prob)
    post_off[data > 0] = 0
    
    post_pois = pois_zero_prob*(1-drop_zero_prob)*(1-off_zero_prob)/((pois_zero_prob*(1-drop_zero_prob) + drop_zero_prob)*(1-off_zero_prob) + off_zero_prob)
    post_pois[data > 0] = 1
    
    print("Computed posterior probabilities.")
    print("Fitting parameters of poisson component...")
    
    # Weighted Count Component
    O = t(matrix(rep(colSums(log(data+epsilon)*(post_pois))/colSums(post_pois),dim(data)[1]),ncol = dim(data)[1]))
    Z = log(data+epsilon) - O

    ALPHA = NULL
    for(j in 1:dim(data)[2]){
      count = data[,j]
      # Ideally intercept is included...
      m1 <- glm(count ~ 0 + W + offset(O[,j]),weights = post_pois[,j],family = poisson(log))
      coefs = m1$coefficients[1:K]
      coefs[is.na(coefs)] = 0
      ALPHA = cbind(ALPHA,coefs)
    }
    leY = W %*% ALPHA + O
    eY = exp(leY)
    
    print("Fitting parameters of drop-out component...")
    
    # Weighted Drop Component
    drop_zero_prob = NULL
    B = NULL
    ALPHA2 = NULL
    for(i in 1:dim(data)[1]){
      success = as.vector(post_drop[i,])
      failure = as.vector(post_pois[i,])
      logit_obj = glm(cbind(success,failure) ~ oW ,family=quasibinomial(logit))
      ALPHA2 = cbind(ALPHA2,logit_obj$coefficients[2:(M+1)])
      B[i] = logit_obj$coefficients[1]
      drop_zero_prob = rbind(drop_zero_prob,logit_obj$fitted.values)
    }
    if(to.plot){
      hist(ALPHA2,breaks = 20)
      print(summary(t(ALPHA2)))
    }
    
    print("Fitting parameters of off component...")
    
    ## Weighted Off Component
    off_zero_prob = NULL
    BI = NULL
    for(j in 1:dim(data)[2]){
      success = as.vector(post_off[,j])
      failure = as.vector(post_drop[,j] + post_pois[,j])
      logit_obj = glm(cbind(success,failure) ~ 1 ,family=quasibinomial(logit))
      BI[j] = logit_obj$coefficients[1]
      off_zero_prob = cbind(off_zero_prob,logit_obj$fitted.values)
    }
    if(to.plot){
      plot(colMeans(data > 0),1-(1/(1+exp(-BI))))
    }
  }
  
  ## Normalization Step 1: Impute
  imputed_data = data 
  imputed_data[data == 0] = (post_drop*eY)[data == 0]
  imputed_data[imputed_data > max(data)] = max(data)
  
  print("Imputed Missing Data.")
  if(!impute.only){
  ## Normalization Step 2: Adjust
  print("Adjusting for Unwanted Factors...")
  O = t(matrix(rep(colMeans(imputed_data),dim(data)[1]),ncol = dim(data)[1]))
  Z = log(imputed_data) - O
  lZ = Z
  for(j in 1:dim(data)[2]){
    count = lZ[,j]
    m1 <- lm(count ~ W,weights = 1-post_off[,j])
    lZ[,j] = post_off[,j]*Z[,j] + (1-post_off[,j])*(m1$coefficients[1] + m1$residuals)
  }
  lZ = lZ + O
  return(t(exp(lZ)))
  }else{
    return(t(imputed_data))
  }
}


