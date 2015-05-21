# Normalization Analogous to http://www.biomedcentral.com/1471-2105/12/480
TechCorrect = function(eSet,gf.vec = NULL, NUM_CELLS_PER_BIN = 10, MAX_NUM_BINS = 10, QTHRESH = .01, MAX_EXP_PCS = 5, MAX_QUAL_PCS = 5, PROP_CUTOFF = .90, corr_method = "spearman", plot.dir = NULL, ignore.zeroes = F,restore.zeroes = T){
  ##-----Correcting for Technical Quality Measures----
  
  # Create plot directory, if necessary
  if (!is.null(plot.dir) && !file.exists(plot.dir)){
    dir.create( plot.dir)
  }
  
  # Gene filter vector
  if(is.null(gf.vec)){
    gf.vec = rep(T,dim(eSet)[1])
  }
  
  # Process Quality Features (Log, Abs-Log, etc...)
  qual = processQf(pData(protocolData(eSet)),rownames(protocolData(eSet)))
  qual = t(na.omit(t(qual)))
  print(summary(qual))
  
  # Initialize Corrected TPM Matrix
  corrected.tpm.matrix = exprs(eSet)
  is.zero = corrected.tpm.matrix  == 0
  
  ##----- Expression PCs (based on high-quality genes)
  pc = prcomp(t(log10(corrected.tpm.matrix[gf.vec,] + 1)), center = TRUE, scale = TRUE)
  
  
  ##----- Quality PCs
  # Selecting relevant features
  qual = qual[,apply(qual,2,sd) > 0]
  
  print("Baseline QC Associations:")
  cors = cor2(pc$x[,1: MAX_EXP_PCS],qual,corr_method = corr_method)
  print (-log10(cors$p))
  
  cors = cor(pc$x[,1: MAX_EXP_PCS],qual,method = corr_method)
  Fish = (1/2)*log((1+cors)/(1-cors))
  z = sqrt((dim(pc$x)[1]-3)/1.06)*Fish
  p = 2*pnorm(-abs(z))
  to.keep = p.adjust(p,method = "BH") < QTHRESH
  to.keep.vec = colSums(matrix(to.keep,nrow =  MAX_EXP_PCS)) > 0
  if (sum(to.keep.vec) == 0){
    print("Nothing to adjust.")
    return(corrected.tpm.matrix)
  }
  keep.quals = qual[,to.keep.vec]
  qpc = prcomp(keep.quals,center = T,scale. = T)
  #qpc=list(x = qual)
  
  correct_names = colnames(qual)[colSums(matrix(to.keep,nrow =  MAX_EXP_PCS)) > 0]
  print(paste("Adjusting for",paste(correct_names,collapse = ","),collapse = " "))
    
  # Technical feature PCA
  max_index = min(min(which(cumsum(qpc$sdev^2) >= PROP_CUTOFF*sum(qpc$sdev^2)), MAX_QUAL_PCS))
  if (!is.null(plot.dir)){
    pdf(paste0(plot.dir,"/qual_cumsum.pdf"))
    plot(cumsum(qpc$sdev^2/sum(qpc$sdev^2)))
    abline(h = PROP_CUTOFF, lty = 2, col = "red")
    dev.off()
  }
  print(paste("Adjusting for",max_index,"axes of technical variation."))
  
  num_x = MAX_EXP_PCS
  num_y = dim(qual)[2]
  
  # Uncorrected correlations between PCs and technical quality features
  if (!is.null(plot.dir)){
    pcvq = cor(pc$x[,1:num_x],qual,method = corr_method)
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
    
  numbins = min(MAX_NUM_BINS,floor((dim(eSet)[2])/NUM_CELLS_PER_BIN))
  print(paste("Num bins set to ",numbins,"."))
  print(paste("Eliminating",max_index*(numbins-1),"degree/s of freedom..."))
  
  
  for (index in max_index:1){
    
    quality = qpc$x[,index]
    
    membership = cut(rank(quality),numbins,labels = F)
    
    all.medians = apply(corrected.tpm.matrix,1,median, na.rm = T) # Gene medians across all samples
    
    # Correction applied to each bin
    for (i in 1:numbins){
      
      bin.samples = (membership == i)
      
      # Bin median for each gene
      gene.medians = apply(corrected.tpm.matrix[,bin.samples],1,median, na.rm = T)
      
      new_values = log10(corrected.tpm.matrix[,bin.samples]+EPSILON) + log10(all.medians + EPSILON) - log10(gene.medians + EPSILON)
      corrected.tpm.matrix[,bin.samples] = 10^new_values - EPSILON
    }
    
    print(paste("Adjusted axis",index))
    
    print("QC Associations before re-zero")
    pc = prcomp(t(log10(corrected.tpm.matrix[gf.vec,] + 1)), center = TRUE, scale = TRUE)
    cors = cor2(pc$x[,1: MAX_EXP_PCS],qual,corr_method = corr_method)
    print (-log10(cors$p[1,]))
    
    # Enforce non-negative TPM
    if(!ignore.zeroes){
      corrected.tpm.matrix[corrected.tpm.matrix < 0] = 0
    }else{
      if(sum(corrected.tpm.matrix < 0,na.rm = T) > 0){
        stop(paste("Generated Negative TPM."))
      }
    }
    
    print("QC Associations after re-zero")
    pc = prcomp(t(log10(corrected.tpm.matrix[gf.vec,] + 1)), center = TRUE, scale = TRUE)
    cors = cor2(pc$x[,1: MAX_EXP_PCS],qual,corr_method = corr_method)
    print (-log10(cors$p[1,]))
  }
  
  if(ignore.zeroes){
    corrected.tpm.matrix[is.na(corrected.tpm.matrix)] = 0
  }
  
  if (!is.null(plot.dir)){
    pdf(paste0(plot.dir,"/adjreads_nz.pdf"))
    plotDensities(log10(corrected.tpm.matrix[gf.vec,] + 1), xlim = c(0,6), ylim = c(0,1), main = "", xlab = "", ylab = "", zero.omit = F)
    dev.off()
    
    pdf(paste0(plot.dir,"/adjcorrs_nz.pdf"))
    pc = prcomp(t(log10(corrected.tpm.matrix[gf.vec & apply(corrected.tpm.matrix,1,sd) > 0,] + 1)), retx = TRUE, center = TRUE, scale = TRUE)
    
    pcvq = cor(pc$x[,1:num_x],qual,method = corr_method)
    Fish = (1/2)*log((1+pcvq)/(1-pcvq))
    z = sqrt((dim(pc$x)[1]-3)/1.06)*Fish
    p = 2*pnorm(-abs(z))
    
    barplot(-log10(t(p)), ylim = c(0,scale),
            main = paste("Z-test of",num_x,"Expression PCs against",num_y,"Quality Feautures"),
            beside = T,xlab = "Expression PCs",
            legend = colnames(qual),args.legend = c( cex = .5))
    dev.off()
    
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
    
    pcvq = cor(pc$x[,1:num_x],qual,method = corr_method)
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