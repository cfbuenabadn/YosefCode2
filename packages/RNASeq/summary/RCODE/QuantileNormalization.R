##----- Extend QN Adjustment Ratios from Normalized Genes to Un-Normalized Genes -----
# norm = quantile normalized values, with NA for all missing normalized data
# raw = un-normalized data

ExtendRat = function(norm,raw){
  
  # Order of raw data
  o = order(raw)
  
  # Ratio Vector (= Normalized Value / Un-Normalized Value)
  r = norm[o]/raw[o]
  is.zero = (raw[o] == 0)
  
  r = r[!is.zero] # Remove Zeroes on the Left
  r[is.infinite(r)] = NA # Infinite ratios are not applicable (will not be used in computing new ratios)
  len = length(r) # Number of ratios for non-zero data                                                            
  
  # Leading NA propagates as -1 (No information gained from this entry)
  if(is.na(r[1])){
    r[1] = -1
  }
  
  # Propagate Ratio to the Right (Define a Left Ratio)
  left.r = r[2:len]
  na.index = which(is.na(left.r))
  while(length(na.index) > 0){
    left.r[na.index] = c(r[1],left.r)[na.index]
    na.index = which(is.na(left.r))
  }
  left.r = c(r[1],left.r)
  left.r[left.r == -1] = NA
  
  # Restore Leading NA
  if(r[1] == -1){
    r[1] = NA
  }
  
  # Tail NA propagates as -1 (N Information Gained From this Entry)
  if(is.na(r[len])){
    r[len] = -1
  }
  
  # Propagate Ratio to the Left (Define a Right Ratio)
  right.r = r[1:(len-1)]
  na.index = which(is.na(right.r))
  while(length(na.index) > 0){
    right.r[na.index] = c(right.r,r[len])[na.index+1]
    na.index = which(is.na(right.r))
  }
  right.r = c(right.r,r[len])
  right.r[right.r == -1] = NA
  
  # Compute ratio from mean of left and right (Add Zeroes Back In!)
  new.r = c(rep(0,sum(is.zero)),rowMeans(cbind(left.r,right.r),na.rm = T))
  
  # Apply ratio
  return((new.r*raw[o])[order(o)])
} 

QuartileNormalization = function(x, gf.vec = NULL, QN.ref.vec = NULL, plot.dir = NULL){
  
  # Create plot directory, if necessary
  if (!is.null(plot.dir) && !file.exists(plot.dir)){
    dir.create(plot.dir)
  }
  
  # Gene filter vector
  if(is.null(gf.vec)){
    gf.vec = rep(T,dim(x)[1])
  }
  
  
  # Upper Quartile Normalization
  quantiles = apply(exprs(tf.sc.eSet)[gf.vec,],2,quantile)
  quartiles = matrix(unlist(quantiles),nrow = 5)[4,]
  if(any(quartiles == 0)){
   stop("Third Quartile Equals Zero! Try more stringent gene/sample filtering.") 
  }
  QRN.matrix = exprs(tf.sc.eSet)
  
  QRN.matrix = t(t(QRN.matrix)/quartiles)
    
  
  if (!is.null(plot.dir)){
    pdf(paste0(plot.dir,"/before_QRN.pdf"))  
    plot(x = NULL, ylim = c(0,1),xlim = c(0,5))
    cols = rainbow(dim(x)[2])
    for(si in 1:dim(x)[2]){
      lines(density(log10(x[gf.vec,][,si]+1)),col = cols[si])
    }
    dev.off()
    
    pdf(paste0(plot.dir,"/before_QRN_all.pdf"))  
    plot(x = NULL, ylim = c(0,1),xlim = c(0,5))
    cols = rainbow(dim(x)[2])
    for(si in 1:dim(x)[2]){
      lines(density(log10(x[,si]+1)),col = cols[si])
    }
    dev.off()
    
    pdf(paste0(plot.dir,"/after_QRN.pdf"))
    plot(x = NULL, ylim = c(0,1),xlim = c(0,5))
    cols = rainbow(dim(x)[2])
    for(si in 1:dim(x)[2]){
      lines(density(log10(QRN.matrix[gf.vec,][,si]+1)),col = cols[si])
    }
    dev.off()
    
    pdf(paste0(plot.dir,"/after_QRN_all.pdf"))
    plot(x = NULL, ylim = c(0,1),xlim = c(0,5))
    cols = rainbow(dim(x)[2])
    for(si in 1:dim(x)[2]){
      lines(density(log10(QRN.matrix[,si]+1)),col = cols[si])
    }
    dev.off()
    
    pdf(paste0(plot.dir,"/bx_before_qrn.pdf"),width = 20)
    boxplot(log10(x[gf.vec,]+1))
    dev.off()
    
    pdf(paste0(plot.dir,"/bx_after_qrn.pdf"),width = 20)
    boxplot(log10(QRN.matrix[gf.vec,]+1))
    dev.off()
    
  }
  
  return(QRN.matrix)
  
}


QuantileNormalization = function(x, gf.vec = NULL, plot.dir = NULL){
  
  # Create plot directory, if necessary
  if (!is.null(plot.dir) && !file.exists(plot.dir)){
    dir.create(plot.dir)
  }
  
  # Gene filter vector
  if(is.null(gf.vec)){
    gf.vec = rep(T,dim(x)[1])
  }
  
  QN.matrix = x
  
  print("Normalizing...")
  
  QN.matrix[gf.vec,] = normalize.quantiles(x[gf.vec,])
  una.QN.matrix = QN.matrix
  # Extending Quantiles and Re-Introducing Zeroes.
  for(si in 1:dim(x)[2]){
    
    z = x[,si]
    z[gf.vec] = QN.matrix[,si][gf.vec]
    z[!gf.vec] = NA
                          
    QN.matrix[,si] = ExtendRat(z,x[,si])
  }
  
  stopifnot(max(QN.matrix[x == 0]) == 0)
  
  if (!is.null(plot.dir)){
    pdf(paste0(plot.dir,"/before_QN.pdf"))  
    plot(x = NULL, ylim = c(0,1),xlim = c(0,5))
    cols = rainbow(dim(x)[2])
    for(si in 1:dim(x)[2]){
      lines(density(log10(x[gf.vec,][,si]+1)),col = cols[si])
    }
    dev.off()
    
    pdf(paste0(plot.dir,"/before_QN_all.pdf"))  
    plot(x = NULL, ylim = c(0,1),xlim = c(0,5))
    cols = rainbow(dim(x)[2])
    for(si in 1:dim(x)[2]){
      lines(density(log10(x[,si]+1)),col = cols[si])
    }
    dev.off()
    
    pdf(paste0(plot.dir,"/after_QN_unzeroed.pdf"))
    plot(x = NULL, ylim = c(0,1),xlim = c(0,5))
    cols = rainbow(dim(x)[2])
    for(si in 1:dim(x)[2]){
      lines(density(log10(una.QN.matrix[gf.vec,][,si]+1)),col = cols[si])
    }
    dev.off()
      
    pdf(paste0(plot.dir,"/after_QN.pdf"))
    plot(x = NULL, ylim = c(0,1),xlim = c(0,5))
    cols = rainbow(dim(x)[2])
    for(si in 1:dim(x)[2]){
      lines(density(log10(QN.matrix[gf.vec,][,si]+1)),col = cols[si])
    }
    dev.off()
    
    pdf(paste0(plot.dir,"/after_QN_all.pdf"))
    plot(x = NULL, ylim = c(0,1),xlim = c(0,5))
    cols = rainbow(dim(x)[2])
    for(si in 1:dim(x)[2]){
      lines(density(log10(QN.matrix[,si]+1)),col = cols[si])
    }
    dev.off()
  }
  
  return(QN.matrix)
  
}
