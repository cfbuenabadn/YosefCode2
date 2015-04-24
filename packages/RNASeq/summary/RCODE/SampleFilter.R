# Data Cleaning Module 2: Filter Samples According to NREADS, RALIGN, Efficiency, and FNR AUC
# Michael Cole, March 2015
# -------------------------

library(diptest)
library(mixtools)

# Function: FNR
# Usage: Fit Logistic Regression Model of FNR against Reference Set of Genes
# Params:
# -------------------------
# eSet = expression set.
# bulk.eSet = expression set. set from which population means are computed
# ref_list = character. path to list of reference gene symbols (one per line)
# FN_thresh = numeric. hard threshold for false negative 
# out.dir = character. path to plot/report directory. No plots if NULL


FNR = function(eSet, bulk.eSet,ref_list, FN_thresh = 0, report_miss = T, out.dir = NULL, plot.name = "FNR.pdf"){
    
    ref.genes = as.character(unlist((read.table(ref_list))))
    
    # Create output directory, if necessary
    if (!is.null(out.dir) && !file.exists(out.dir)){
      dir.create(out.dir)
    }
    
    common_symbols = gsub("_variant[[:digit:]]*","",featureData(eSet)$Gene_Symbol)
    is.shared = ref.genes %in% common_symbols
    
    # Report missing reference genes
    num_missing = sum(!is.shared)
    if (report_miss && (num_missing  > 0)){
      if (!is.null(out.dir)){
        warning(paste("Warning:",num_missing,"reference transcripts missing from expression set: writing to missing_reference.txt"))
        write.table(ref.genes[!is.shared],file = paste0(out.dir,"/missing_reference.txt"),sep = "\t",quote = F,row.names = F, col.names = F)
      }else {
        warning(paste("Warning:",num_missing,"reference transcripts missing from expression set."))
      }
    }
    
    is.Reference = common_symbols %in% ref.genes[is.shared]
    
    stopifnot(!any(rownames(eSet) != rownames(bulk.eSet))) # Stop if bulk set has different rows (transcripts)
    ref.bulk.expr.vec = log10(rowMeans(exprs(bulk.eSet)[is.Reference,])+1)
    ref.drop.out.matrix = (exprs(eSet)[is.Reference,] <= FN_thresh)
    
    # Logistic Regression Model of FNR
    ref.glms = list()
    for (si in 1:dim(ref.drop.out.matrix)[2]){
      fit = suppressWarnings(glm(cbind(ref.drop.out.matrix[,si],1 - ref.drop.out.matrix[,si]) ~ ref.bulk.expr.vec,family=binomial(logit)))
      if(fit$converged){
        ref.glms[[si]] = fit$coefficients
      } else {
        ref.glms[[si]] = NA
      }
    }
    
    # Plot FNR Curves
    if (!is.null(out.dir)){
      pdf(paste0(out.dir,"/",plot.name))
      plot(x = NULL, xlim = c(0,6), ylim = c(0,1), ylab = "1-FNR",xlab = "log10(Mean TPM+1) in Bulk", main = "FNR in Reference Genes")
      for (si in 1:dim(ref.drop.out.matrix)[2]){
        if(!any(is.na(ref.glms[[si]]))){
          x = seq(0,6,.01)
          y = colSums(ref.glms[[si]] * rbind(rep(1,length(x)),x))
          z = 1 - exp(y) /(1 + exp(y))
          lines(x,z)
        }
      }
      dev.off()
    }
    
    # Return Fits
    return(ref.glms)
    
}

FNRw = function(eSet, bulk.eSet, gf.vec, FN_thresh, housekeeping_list){
  stopifnot(!any(rownames(eSet) != rownames(bulk.eSet))) # Stop if bulk set has different rows (transcripts)
  ref.glms = FNR(eSet[gf.vec,], bulk.eSet[gf.vec,],housekeeping_list, FN_thresh, report_miss = F, out.dir = NULL)
  
  
  drop.out.matrix = (exprs(eSet) <= FN_thresh)
  bulk.expr.vec = log10(rowMeans(exprs(bulk.eSet))+1)
  weights = 1-drop.out.matrix
  
  for (si in 1:dim(weights)[2]){
    x = bulk.expr.vec[drop.out.matrix[,si]]
    y = colSums(ref.glms[[si]] * rbind(rep(1,length(x)),x))
    z = exp(y) /(1 + exp(y))
    weights[,si][drop.out.matrix[,si]] = z
  }
  rownames(weights) = featureData(eSet)$Gene
  
  return(weights)
  
}

# ref_list = character. path to list of reference gene symbols (one per line)

SampleFilter = function(eSet, gene.filter.vec = NULL, housekeeping_list, mixture = T, verbose = T, plot.dir = NULL){
  
  # Create plot directory, if necessary
  if (!is.null(plot.dir) && !file.exists(plot.dir)){
    dir.create(plot.dir)
  }
  
  qual.names = c("NREADS","RALIGN")
  
  ##----- Hard Cutoffs
  MIN_NREADS = 25000
  MIN_RALIGN = 15 # Used to be 20
  MIN_EFF = 0.2
  
  ##----- Sufficient Cutoffs
  SUFF_RALIGN = 65
  SUFF_EFF = 0.8
  
  ##----- Deviation
  Z_CUTOFF = 1.645
  
  ##----- AUC Params
  MAX_LOG = 6 
  
  ##----- Calculate Final Cutoffs
  
  # Number of reads (Transformed)
  logr = unlist(log10(pData(protocolData(eSet)[,qual.names[1]]) + 1))
  
  # Simple thresholds
  LOGR_CUTOFF = max(median(logr) - Z_CUTOFF*mad(logr), log10(MIN_NREADS + 1))
  LOGR_CUTOFF = max(mean(logr) - Z_CUTOFF*sd(logr), LOGR_CUTOFF)
  
  # Mixture model threshold
  if(mixture){    
    
    is.multimodal = dip.test(logr)$p.value < .01
    
    if(is.multimodal){
      print("NREADS is multimodal -> Applying normalmixEM...")
      mixmdl = normalmixEM(logr,k=2)
      component = which(mixmdl$mu %in% max(mixmdl$mu))
      LOGR_CUTOFF = max(mixmdl$mu[component] - Z_CUTOFF*mixmdl$sigma[component], LOGR_CUTOFF)
    }
  }
  
  # Ratio of reads aligned
  ralign = unlist(pData(protocolData(eSet)[,qual.names[2]]))
  
  # Simple thresholds
  RALIGN_CUTOFF = max(median(ralign) - Z_CUTOFF*mad(ralign), MIN_RALIGN)
  RALIGN_CUTOFF = max(mean(ralign) - Z_CUTOFF*sd(ralign), RALIGN_CUTOFF)
  
  # Mixture model threshold
  if(mixture){    
    
    is.multimodal = dip.test(ralign)$p.value < .01
    
    if(is.multimodal){
      print("RALIGN is multimodal -> Applying normalmixEM...")
      mixmdl = normalmixEM(ralign,k=2)
      component = which(mixmdl$mu %in% max(mixmdl$mu))
      RALIGN_CUTOFF = max(mixmdl$mu[component] - Z_CUTOFF*mixmdl$sigma[component], RALIGN_CUTOFF)
    }
  }
  
  RALIGN_CUTOFF = min(RALIGN_CUTOFF,SUFF_RALIGN)
    
  # Efficiency: Fraction of filtered genes expressed higher than minimum level (typically zero)
  if (is.null(gene.filter.vec)){
    gene.filter.vec = rep(T,dim(eSet)[1])
  }
  efficiency = colMeans(exprs(eSet)[gene.filter.vec,] > min(exprs(eSet)))
  
  # Simple thresholds
  EFF_CUTOFF = max(median(efficiency) - Z_CUTOFF*mad(efficiency), MIN_EFF)
  EFF_CUTOFF = max(mean(efficiency) - Z_CUTOFF*sd(efficiency), EFF_CUTOFF)
  
  # Mixture model threshold
  if(mixture){    
    
    is.multimodal = dip.test(efficiency)$p.value < .01
    
    if(is.multimodal){
      print("Efficiency is multimodal -> Applying normalmixEM...")
      mixmdl = normalmixEM(efficiency,k=2)
      component = which(mixmdl$mu %in% max(mixmdl$mu))
      EFF_CUTOFF = max(mixmdl$mu[component] - Z_CUTOFF*mixmdl$sigma[component], EFF_CUTOFF)
    }
  }
  
  EFF_CUTOFF = min(EFF_CUTOFF,SUFF_EFF)
  
  # FNR AUC Filter
  
  ref.glms = FNR(eSet = eSet[gene.filter.vec,],bulk.eSet = eSet[gene.filter.vec,],housekeeping_list,out.dir = plot.dir,plot.name = "fnr_before.pdf")

  # Compute AUC  
  AUC = NULL
  for (si in 1:dim(sc.eSet)[2]){
    if(!any(is.na(ref.glms[[si]]))){
      y = sum(ref.glms[[si]] * c(1,MAX_LOG))
      AUC[si] = MAX_LOG - log(exp(y) + 1)/ref.glms[[si]][2] + log(exp(ref.glms[[si]][1]) + 1)/ref.glms[[si]][2] 
    } else {
      AUC[si] = NA
    }
  }
  
  is.good_AUC = !is.na(AUC)
  is.FNR_nonconverge = is.na(AUC)
  is.multimodal = dip.test(AUC)$p.value < .01
  
  if(is.multimodal){
    print("AUC is multimodal -> Applying normalmixEM...")
    mixmdl = normalmixEM(na.omit(AUC),k = 2)
    
    # Selecting higher AUC component based on maximum posterior discriminator
    component = which(mixmdl$mu %in% max(mixmdl$mu))
    is.good_AUC[!is.na(AUC)] = mixmdl$posterior[,component] > 0.5
  }
  
  is.Low.Reads = logr < LOGR_CUTOFF
  is.Low.Read.Rat = ralign < RALIGN_CUTOFF
  is.Low.Efficiency = efficiency < EFF_CUTOFF
  is.Low.Technical.Quality = is.Low.Reads | is.Low.Read.Rat | is.Low.Efficiency | !is.good_AUC
  
  sf.eSet = eSet[,!is.Low.Technical.Quality]
  
  if (!is.null(plot.dir)){
    
    report = cbind(colnames(eSet),is.Low.Reads,is.Low.Read.Rat,is.Low.Efficiency,is.FNR_nonconverge,!is.good_AUC)    
    report = report[is.Low.Technical.Quality,]
    colnames(report) = c("Filtered_Sample","Low_Reads","Low_Aligned_Ratio","Low_Efficiency","FNR_Not_Converged","Low/Failed_FNR_AUC")
    write.table(report,file = paste0(plot.dir,"/filter_report.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
    
    pdf(paste0(plot.dir,"/nreads_before.pdf"))
    hist(logr, main = "NREADS Filter", xlab = "log10(NREADS+1)")
    abline(v = LOGR_CUTOFF, col = "red", lty = 2)
    print(paste("NREADS Cutoff:",10^LOGR_CUTOFF - 1))
    dev.off()
    
    pdf(paste0(plot.dir,"/ralign_before.pdf"))
    hist(ralign, main = "RALIGN Filter", xlab = "RALIGN")
    abline(v = RALIGN_CUTOFF, col = "red", lty = 2)
    print(paste("RALIGN Cutoff:",RALIGN_CUTOFF))
    dev.off()
    
    pdf(paste0(plot.dir,"/effic_before.pdf"))
    hist(efficiency, main = "Efficiency Filter", xlab = "Efficiency")
    abline(v = EFF_CUTOFF, col = "red", lty = 2)
    print(paste("Efficiency:",EFF_CUTOFF))
    dev.off()
    
    pdf(paste0(plot.dir,"/nreads_after.pdf"))
    hist(logr[!is.Low.Technical.Quality], main = "NREADS Filter", xlab = "log10(NREADS+1)")
    abline(v = LOGR_CUTOFF, col = "red", lty = 2)
    dev.off()
    
    pdf(paste0(plot.dir,"/ralign_after.pdf"))
    hist(ralign[!is.Low.Technical.Quality], main = "RALIGN Filter", xlab = "RALIGN")
    abline(v = RALIGN_CUTOFF, col = "red", lty = 2)
    dev.off()
    
    pdf(paste0(plot.dir,"/effic_after.pdf"))
    hist(efficiency[!is.Low.Technical.Quality], main = "Efficiency Filter", xlab = "Efficiency")
    abline(v = EFF_CUTOFF, col = "red", lty = 2)
    dev.off()
  }
  
  FNR(eSet = sf.eSet[gene.filter.vec,],bulk.eSet = sf.eSet[gene.filter.vec,],housekeeping_list,report_miss = F,out.dir = plot.dir,plot.name = "fnr_after.pdf")
  
  # Report
  if(verbose){
    print(paste(dim(sf.eSet)[2],"samples out of",dim(eSet)[2],"passing technical quality filter."))
  }
  
  return(sf.eSet)
}


