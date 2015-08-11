require(gplots)

##----- Function for processing raw quality data frame -----
# Transforms quality features to redu
processQf = function(raw.quality.data, row.names, to.log = c("NREADS", "NALIGNED"), to.abs.log = c("MEDIAN_5PRIME_TO_3PRIME_BIAS")){
    
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


# Cell filtering According to QC Summary Score
TechFilter = function(eSet,gf.vec = NULL, Z_CUTOFF = 2.3, mixture = T,PROP_CUTOFF = .90, plot.dir = NULL, plot.prefix = NULL, MAX_QUAL_PCS = 5, MAX_EXP_PCS = 5, MIN_VAR = .7, good.metrics = c("NREADS","NALIGNED","RALIGN","PCT_CODING_BASES") , force.metrics = NULL, DIP_THRESH = .01){
  
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
  #turn Infs to NA so that they will be filtered in the next line. Must to like that because is.infinite does not accept a list (while data.frame is a list)
  qual <- data.frame(lapply(qual, function(x) replace(x, is.infinite(x),NA)), row.names=rownames(qual))
  qual = t(na.omit(t(qual)))
  tpm.matrix = exprs(eSet)
 
  ##----- Expression PCs (based on high-quality genes)
  pc = prcomp(t(log10(tpm.matrix[gf.vec,] + 1)), center = TRUE, scale = TRUE)
  
  
  ##----- Quality PCs
  # Selecting relevant features
  qual = qual[,apply(qual,2,sd) > 0]
  cors = cor(pc$x[,1: MAX_EXP_PCS],qual,method = "spearman")
  Fish = (1/2)*log((1+cors)/(1-cors))
  z = sqrt((dim(pc$x)[1]-3)/1.06)*Fish
  p = 2*pnorm(-abs(z))
  to.keep = p.adjust(p,method = "BH") < .01
  to.keep.vec = colSums(matrix(to.keep,nrow =  MAX_EXP_PCS)) > 0
  print(paste("Selected:",paste(colnames(qual)[to.keep.vec],collapse = " ")))
  
  # Introduce forced metrics
  if(!is.null(force.metrics)){
    to.keep.vec =  to.keep.vec | (colnames(qual) %in% force.metrics)
    print(paste("Forced Selection:",paste(colnames(qual)[to.keep.vec],collapse = " ")))
  }
  
  if (sum(to.keep.vec) == 0){
    print("Nothing to filter.")
    return(eSet)
  }
  keep.quals = qual[,to.keep.vec]
  qpc = prcomp(keep.quals,center = T,scale. = T)
  
  pdf(paste0(plot.dir,"/",plot.prefix,"qual_cumsum.pdf"))
  csum = cumsum((qpc$sdev^2)/sum(qpc$sdev^2))
  plot(csum, main = "Cumulative Quality PC Variance", ylab = "Fraction of Total Variance")
  abline(h = MIN_VAR, lty = 2, col = "red")
  dev.off()
  
  MAX_QUAL_PCS = min(MAX_QUAL_PCS,dim(qpc$x)[2])
  for (i in 1:MAX_QUAL_PCS){
    pdf(paste0(plot.dir,paste0("/",plot.prefix,"qc_pc",i,".pdf")))
    par(mfrow = c(2,1))
    hist(qpc$x[,i],breaks = 20, main = paste0("Distribution of Quality PC ",i), xlab = paste0("Qual PC",i))
    barplot(abs(qpc$rotation[,i]),col = c("red","green")[1 + (qpc$rotation[,i] > 0)], cex.names = .25,horiz = T, las=1, main = "Loadings")
    dev.off()
    if(csum[i] < MIN_VAR){break()}
  }

  
  pdf(paste0(plot.dir,"/",plot.prefix,"qual_corr_heatmap.pdf"))
  heatmap.2(cor(keep.quals),key.title = "",key.xlab = "Spearman Corr.",density.info = 'none',trace = 'none',margin = c(20,20), cexRow = .7, cexCol = .7)
  dev.off()
  
  # Only perfom filtering when Z_CUTOFF is not null
  if (!is.null(Z_CUTOFF)){
    
    # Initializing sample removal vector
    to.remove = rep(F,dim(eSet)[2])
    
    # Check if "good" metrics have been selected -> if filtering is signed
    is.signed = any(colnames(keep.quals) %in% good.metrics)
    if(is.signed){
      print("Signed Filtering.")
    }else{
      print("Unsigned Filtering.")
    }
    
    # Loop over quality PCs until we've covered MIN_VAR of the quality variance
    for ( i in 1:MAX_QUAL_PCS){
      
      # Quality Score
      qscore = qpc$x[,i] 
      
      # Signed Quality Score
      if (is.signed){
        qscore = qscore * median(sign(qpc$rotation[,i][colnames(keep.quals) %in% good.metrics]))
      }
      
      # Simple thresholds
      CUTOFF = median(qscore) - Z_CUTOFF*mad(qscore)
      CUTOFF = max(mean(qscore) - Z_CUTOFF*sd(qscore), CUTOFF)
      
      
      # Reverse thresholds
      if(!is.signed){
        RCUTOFF = median(qscore) + Z_CUTOFF*mad(qscore)
        RCUTOFF = min(mean(qscore) + Z_CUTOFF*sd(qscore), RCUTOFF)
      }
      
      # Mixture model threshold (Signed Only!)
      if(mixture && is.signed){    
        
        is.multimodal = dip.test(qscore)$p.value < DIP_THRESH    
        if(is.multimodal){
          print(paste0("Multimodality detected in Quality Score ",i))
          mixmdl = normalmixEM(qscore,k=2)
          plot(mixmdl,2)
          component = which(mixmdl$mu %in% max(mixmdl$mu))
          CUTOFF = max(mixmdl$mu[component] - Z_CUTOFF*mixmdl$sigma[component], CUTOFF)
          
        }
      }
      
      pdf(paste0(plot.dir,paste0("/",plot.prefix,"qc_filter_pc",i,".pdf")),width = 20)
      par(mfrow = c(1,2))
      hist(qscore,breaks = 20,main = paste0("Distribution of Quality Score ",i),  xlab = paste0("Qual Score ",i))
      abline(v = CUTOFF, lty = 2, col = "red")
      if(!is.signed){
        abline(v = RCUTOFF, lty = 2, col = "red")
      }
      
      if(csum[i] < MIN_VAR){
        if(is.signed){
          to.remove = to.remove | (qscore < CUTOFF)
          hist(qpc$x[,i][!(qscore < CUTOFF)],breaks = 20, main = paste0("Thresh = ",signif(CUTOFF,3),", Rm = ",sum(qscore < CUTOFF),", Tot Rm = ",sum(to.remove) ),xlab = paste0("Qual Score ",i))
        }else{
          to.remove = to.remove | ((qscore < CUTOFF) | (qscore > RCUTOFF))
          hist(qpc$x[,i][!((qscore < CUTOFF) | (qscore > RCUTOFF))],breaks = 20, main = paste0("Threshs = ",signif(CUTOFF,3),"/",signif(RCUTOFF,3),", Rm = ",sum((qscore < CUTOFF) | (qscore > RCUTOFF)),", Tot Rm = ",sum(to.remove) ),xlab = paste0("Qual Score ",i))
        }
        dev.off()
        break()
      }
      else{
        if(is.signed){
          to.remove = to.remove | (qscore < CUTOFF)
          hist(qpc$x[,i][!(qscore < CUTOFF)],breaks = 20, main = paste0("Thresh = ",signif(CUTOFF,3),", Rm = ",sum(qscore < CUTOFF) ),xlab = paste0("Qual Score ",i))
        }else{
          to.remove = to.remove | ((qscore < CUTOFF) | (qscore > RCUTOFF))
          hist(qpc$x[,i][!((qscore < CUTOFF) | (qscore > RCUTOFF))],breaks = 20, main = paste0("Threshs = ",signif(CUTOFF,3),"/",signif(RCUTOFF,3),", Rm = ",sum((qscore < CUTOFF) | (qscore > RCUTOFF)) ),xlab = paste0("Qual Score ",i))
        }
        dev.off()
      }
           
      
    }
    return(eSet[,!to.remove])
  }else{
    return(eSet)
  }
}
