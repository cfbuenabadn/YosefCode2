CalcSig = function(x,w, table, scale = T){
  if (scale){
    x = t(scale(t(x)))
  }
  # Read Sig Table
  sig = read.table(table)
  common_symbols = gsub("_variant[[:digit:]]*","",rownames(x))
  
  sig = sig[sig$V3 %in% common_symbols,]
  sig.names = unique(sig$V1)
  sig.values = matrix(0,nrow = length(sig.names),ncol = dim(x)[2] )
  rownames(sig.values) = sig.names
  
  # For each Signature
  for(i in 1:length(sig.names)){
    sig.indices = which(sig$V1 == sig.names[i])
    sig.signs = c(0,-1,1)[sig$V2[sig.indices]]
    sig.genes = as.character(sig$V3[sig.indices])
    
    # Positive Sum
    is.pos = sig.signs > 0
    is.pos.mapped = common_symbols %in% sig.genes[is.pos]
    num.pos = sum(is.pos.mapped)
    if(num.pos > 1){
      pos.sig.values = colSums(x[is.pos.mapped,]*w[is.pos.mapped,])/colSums(w[is.pos.mapped,])
    }else if(num.pos == 1){
      pos.sig.values = x[is.pos.mapped,]
    }else{
      pos.sig.values = 0 
    }
    
    # Negative Sum
    is.neg = sig.signs < 0
    is.neg.mapped = common_symbols %in% sig.genes[is.neg]
    num.neg = sum(is.neg.mapped)
    if(num.neg > 1){
      neg.sig.values = colSums(x[is.neg.mapped,]*w[is.neg.mapped,])/colSums(w[is.neg.mapped,])
    }else if(num.neg == 1){
      neg.sig.values = x[is.neg.mapped,]
    }else{
      neg.sig.values = 0 
    }
    
    sig.values[i,] = pos.sig.values - neg.sig.values
    
    if (i %% 100 == 0){
    print(paste(i,"/",length(sig.names),"signatures computed"))
    }
  }
  rownames(sig.values) = sig.names
  colnames(sig.values) = colnames(x)
  print("All signatures computed!")
  return(sig.values)
}

GetSigSet = function(table, name){
  sig = read.table(table)
  sig.indices = which(sig$V1 == name)
  sig.signs = c(0,-1,1)[sig$V2[sig.indices]]
  sig.genes = as.character(sig$V3[sig.indices])
  return(cbind(sig.genes,sig.signs))
}


