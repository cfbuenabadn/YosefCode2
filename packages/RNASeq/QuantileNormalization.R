ExtendRat = function(norm,raw){
  o = order(raw)
  
  # Ratio vector
  r = norm[o]/raw[o]
  is.zero = (raw[o] == 0)
  r = r[!is.zero]
  r[is.infinite(r)] = NA
  len = length(r)                                                                 
  
  # Leading NA propagates as -1
  if(is.na(r[1])){
    r[1] = -1
  }
  
  # Left Ratio
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
  
  # Tail NA propagates as -1
  if(is.na(r[len])){
    r[len] = -1
  }
  
  # Right Ratio
  right.r = r[1:(len-1)]
  na.index = which(is.na(right.r))
  while(length(na.index) > 0){
    right.r[na.index] = c(right.r,r[len])[na.index+1]
    na.index = which(is.na(right.r))
  }
  right.r = c(right.r,r[len])
  right.r[right.r == -1] = NA
  
  # Compute ratio from mean of left and right
  new.r = c(rep(0,sum(is.zero)),rowMeans(cbind(left.r,right.r),na.rm = T))
  
  # Apply ratio
  return((new.r*raw[o])[order(o)])
} 


QuantileNormalization = function(x, gf.vec = NULL, QN.ref.vec = NULL, plot.dir = NULL){
  
  # Create plot directory, if necessary
  if (!is.null(plot.dir) && !file.exists(plot.dir)){
    dir.create(plot.dir)
  }
  
  # Gene filter vector
  if(is.null(gf.vec)){
    gf.vec = rep(T,dim(x)[1])
  }
  
  QN.matrix = x
  
  # Reference vector for QN
  if (is.null(QN.ref.vec)){
    QN.ref.matrix = x[gf.vec,] 
    for(si in 1:dim(x)[2]){
      QN.ref.matrix[,si] = sort(QN.ref.matrix[,si])
    }
    QN.ref.vec = rowMeans(QN.ref.matrix)
  }
  
  print("Normalizing...")
  
  # Normalization
  for(si in 1:dim(x)[2]){
    
    # Expression Vector of Subject
    samp.m.vec = x[,si][gf.vec]
    
    # Unadjusted QN
    comp.m.vec = sort(QN.ref.vec)
    comp.m.vec[sort(samp.m.vec) == 0] = 0
    comp.m.vec = comp.m.vec[order(order(samp.m.vec))]
    
    z = x[,si]
    z[gf.vec] = comp.m.vec
    z[!gf.vec] = NA
    
    QN.matrix[,si] = ExtendRat(z,x[,si])
  }
  
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