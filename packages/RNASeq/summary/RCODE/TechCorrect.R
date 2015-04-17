##----- Function for processing raw quality data frame -----
processQf = function(raw.quality.data, row.names){
  
  to.log = c("NREADS", "NALIGNED")
  to.abs.log = c("MEDIAN_5PRIME_TO_3PRIME_BIAS")
  
  quality.features = raw.quality.data
  rownames(quality.features) = row.names
  
  for (i in 1:dim(quality.features)[2]){
    if(!is.numeric(quality.features[,i])){
      quality.features[,i] = as.numeric(as.character(quality.features[,i]))
    }
  }
  
  # Special Transformations
  quality.features[,to.log] = log10(quality.features[,to.log])
  quality.features[,to.abs.log] = exp(abs(log(quality.features[,to.abs.log])))
  
  return(quality.features)
}

# Normalization Analogous to http://www.biomedcentral.com/1471-2105/12/480
TechCorrect = function(eSet,gf.vec = NULL, maxnumbins = 5, Z_CUTOFF = 1, PROP_CUTOFF = .90, plot.dir = NULL, ignore.zeroes = F,restore.zeroes = T){
    
  ##-----Correcting for Technical Quality Measures----
  
  # Create plot directory, if necessary
  if (!is.null(plot.dir) && !file.exists(plot.dir)){
    dir.create( plot.dir)
  }
  
  # Gene filter vector
  if(is.null(gf.vec)){
    gf.vec = rep(T,dim(eSet)[1])
  }
  
  qual = processQf(pData(protocolData(eSet)),rownames(protocolData(eSet)))
  qual = t(na.omit(t(qual)))
  corrected.tpm.matrix = exprs(eSet)
  is.zero = corrected.tpm.matrix  == 0
  
  # Bin Params
  binsize = ceiling(dim(corrected.tpm.matrix)[2]/maxnumbins)
  numbins = ceiling(dim(corrected.tpm.matrix)[2]/binsize)
  while( dim(corrected.tpm.matrix)[2]-binsize*(numbins-1) < binsize/2){
    maxnumbins = maxnumbins-1
    if(maxnumbins == 0){
      stop("Can't space binning. Increase sample size ;)")
      return(NULL)
    }
    binsize = ceiling(dim(corrected.tpm.matrix)[2]/maxnumbins)
    numbins = ceiling(dim(corrected.tpm.matrix)[2]/binsize)
  }
  
  print(paste("Number of bins set to",numbins))
  
  ##----- Expression PCs (based on high-quality genes)
  pc = prcomp(t(log10(corrected.tpm.matrix[gf.vec,] + 1)), center = TRUE, scale = TRUE)
  
  
  ##----- Quality PCs
  # Selecting relevant features
  
  cors = cor(pc$x,qual,method = "spearman")
  scores = rev(sort(colSums(pc$sdev^2/sum(pc$sdev^2) * cors^2)))
  to.correct = names(scores)[scores > median(scores) + Z_CUTOFF*mad(scores)]
  print(paste("Adjusting for",paste(to.correct,collapse = ","),collapse = " "))
  
  # Show features selected
  if (!is.null(plot.dir)){
    pdf(paste0(plot.dir,"/feature_scores.pdf"))
    plot(x = NULL, ylim = c(0,max(scores)),xlim = c(0,length(scores)+1), ylab = "Correlation Score", main = "Technical Feature Selection")
    text(sort(scores),labels = names(scores)[order(scores)], col = c("black","red")[1+as.numeric(sort(scores) > median(scores) + Z_CUTOFF*mad(scores))])
    abline(h = median(scores) + Z_CUTOFF*mad(scores), lty = 2, col = "red")
    dev.off()
  }
  
  # Technical feature PCA
  qpc = prcomp(qual[,to.correct],center = T,scale. = T)
  max_index = min(which(cumsum(qpc$sdev^2) >= PROP_CUTOFF*sum(qpc$sdev^2)))
  if (!is.null(plot.dir)){
    pdf(paste0(plot.dir,"/tech_cumcont.pdf"))
    plot(cumsum(qpc$sdev^2/sum(qpc$sdev^2)))
    abline(h = PROP_CUTOFF, lty = 2, col = "red")
    dev.off()
  }
  print(paste("Adjusting for",max_index,"axes of technical variation."))
  print(paste("Eliminating",max_index*(numbins-1),"degree/s of freedom..."))
  
  MAX_NUMX = 5
  num_x = min(max(which(pc$sdev^2 >= sum(pc$sdev^2)/length(pc$sdev))), MAX_NUMX)
  num_y = dim(qual)[2]
  # Uncorrected correlations between PCs and technical quality features
  if (!is.null(plot.dir)){
    pcvq = cor(pc$x[,1:num_x],qual,method = "spearman")
    Fish = (1/2)*log((1+pcvq)/(1-pcvq))
    z = sqrt((dim(pc$x)[1]-3)/1.06)*Fish
    p = 2*pnorm(-abs(z))
    scale = max(-log10(p))
    
    pdf(paste0(plot.dir,"/rawcorrs.pdf"))
    barplot(-log10(t(p)), ylim = c(0,scale),
            main = paste("Z-test of",num_x,"Expression PCs against",num_y,"Quality Feautures"),
            beside = T,xlab = "Expression PCs",
            legend = colnames(qual),args.legend = c( cex = .5))
    dev.off()
  }
  
  # Correct for each significant technical PC
  
  if(!ignore.zeroes){
    EPSILON = 1
  }else{
    EPSILON = 0
    corrected.tpm.matrix[corrected.tpm.matrix == 0] = NA
  }
  
  for (index in max_index:1){
    quality = qpc$x[,index]
    names(quality) = rownames(qpc$x)
    all.means = apply(corrected.tpm.matrix,1,mean, na.rm = T) # Gene means across all samples
    
    # Correction applied to each bin
    for (i in 1:numbins){
      bin.samples = names(sort(quality)[(1 + (i-1)*binsize):min(binsize + (i-1)*binsize,dim(corrected.tpm.matrix)[2])])
      # Bin mean for each gene
      gene.means = apply(corrected.tpm.matrix[,bin.samples],1,mean, na.rm = T)
      new_values = log10(corrected.tpm.matrix[,bin.samples]+EPSILON) + log10(all.means + EPSILON) - log10(gene.means + EPSILON)
      corrected.tpm.matrix[,bin.samples] = 10^new_values - EPSILON
    }
    print(paste("Adjusted axis",index))
    
    # Enforce non-negative TPM
    if(!ignore.zeroes){
      corrected.tpm.matrix[corrected.tpm.matrix < 0] = 0
    }else{
      if(sum(corrected.tpm.matrix < 0,na.rm = T) > 0){
        stop(paste("Generated Negative TPM."))
      }
    }
  }
  
  if(ignore.zeroes){
    corrected.tpm.matrix[is.na(corrected.tpm.matrix)] = 0
  }
  if(restore.zeroes){
    corrected.tpm.matrix[is.zero] = 0
  }
  
  #Expression PCs vs Quality Features
  if (!is.null(plot.dir)){
    pdf(paste0(plot.dir,"/adjreads.pdf"))
    plotDensities(log10(corrected.tpm.matrix[gf.vec,] + 1), xlim = c(0,6), ylim = c(0,1), main = "", xlab = "", ylab = "", zero.omit = F)
    dev.off()
    
    pdf(paste0(plot.dir,"/adjcorrs.pdf"))
    pc = prcomp(t(log10(corrected.tpm.matrix[gf.vec & apply(corrected.tpm.matrix,1,sd) > 0,] + 1)), retx = TRUE, center = TRUE, scale = TRUE)
    
    pcvq = cor(pc$x[,1:num_x],qual,method = "spearman")
    Fish = (1/2)*log((1+pcvq)/(1-pcvq))
    z = sqrt((dim(pc$x)[1]-3)/1.06)*Fish
    p = 2*pnorm(-abs(z))
    
    barplot(-log10(t(p)), ylim = c(0,scale),
            main = paste("Z-test of",num_x,"Expression PCs against",num_y,"Quality Feautures"),
            beside = T,xlab = "Expression PCs",
            legend = colnames(qual),args.legend = c( cex = .5))
    dev.off()
    
    pdf(paste0(plot.dir,"/adjrat.pdf"))
    plot(rowMeans(exprs(eSet)[gf.vec,] == 0),log10(rowMeans(corrected.tpm.matrix[gf.vec,]/exprs(eSet)[gf.vec,],na.rm = T)),xlab = "FNR",ylab = "log10(Mean-Adjustment)")
    dev.off()
  }
  
  return(corrected.tpm.matrix)
}