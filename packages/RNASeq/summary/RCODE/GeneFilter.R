# Data Cleaning Module 1: Filter Low-TPM Genes + Visualization
# Michael Cole, March 2015
# -------------------------

# Function: plotDensities
# Usage: Plot read densities
# Params:
# -------------------------
# x = numeric matrix. Missing values are NOT allowed
# xlim, ylim, main, xlab, ylab = plot parameters
# zero.omit = logical. If true, will omit zero reads from plot

plotDensities = function(x, xlim, ylim, main, xlab, ylab, zero.omit = F){
  if(zero.omit){
    x[x == 0] = NA
  }
  plot(x = NULL,
       xlim = xlim, ylim = ylim, 
       ylab = ylab, 
       xlab = xlab, 
       main = main)
  
  # Add one density per sample (Colors of the rainbow)
  cols = rainbow(dim(x)[2])
  density.list = list()
  for(i in 1:dim(x)[2]){
    density.list[[i]] = density(x[,i], na.rm = T)
    lines(density.list[[i]], col = cols[i])
  }
}

# Function: GeneFilter
# Usage: Produce Logical Gene Filter Vector and Visualize Cuts
# Params:
# -------------------------
# eSet = expression set
# count.cutoff = numeric. low-count threshold
# prop.failed = numeric. gene is filtered if propfailed or more samples fall below the low-count threshold
# verbose = logical. If true, will print message describing filtering step
# plot.dir = character. path to plot directory. No plots if NULL

GeneFilter = function(eSet, count.cutoff, prop.failed, verbose = F, plot.dir = NULL, plot.prefix = NULL, is_abs = F){
  
  # Create plot directory, if necessary
  if (!is.null(plot.dir) && !file.exists(plot.dir)){
    dir.create(plot.dir)
  }
  
  # Plot Raw Read distribution
  if (!is.null(plot.dir)){
    log.tpm.matrix = log10(exprs(eSet) + 1)
    pdf(paste0(plot.dir,"/",plot.prefix,"_rawreads.pdf"))
    plotDensities( x = log.tpm.matrix,
                   xlim = c(0,6), ylim = c(0,1),
                   main = "Read Distribution - Raw",
                   xlab = "log10(TPM + 1)",
                   ylab = "Read Kernel Density")
    abline(v = log10(count.cutoff+1), lty = 2, lwd = 2.5, col = "Red")
    dev.off()
  }
  
  # Implement cut
  if(!is_abs) {
    #relative (original) filter
    is.Cut.Gene = rowMeans(exprs(eSet) < count.cutoff) >= prop.failed
  } else {
    #absolute filter
    is.Cut.Gene = rowSums(exprs(eSet) < count.cutoff) >= prop.failed
  }
    
  #quick&dirty fix to cut by absolute values
  #is.Cut.Gene = !(rowSums(exprs(eSet) >= count.cutoff) >= 20) #has at least 20 cells in which expressed
  genefilter.eSet = eSet[!is.Cut.Gene,]
  
  # Plot Cut Read distribution
  if (!is.null(plot.dir)){
    log.cut.tpm.matrix = log10(exprs(genefilter.eSet) + 1)
    pdf(paste0(plot.dir,"/",plot.prefix,"_cutreads.pdf"))
    plotDensities( x = log.cut.tpm.matrix,
                   xlim = c(0,6), ylim = c(0,1),
                   main = "Read Distribution - Cut",
                   xlab = "log10(TPM + 1)",
                   ylab = "Read Kernel Density")
    abline(v = log10(count.cutoff+1), lty = 2, lwd = 2.5, col = "Red")
    dev.off()
  }
  
  if(verbose){
    print(paste0("Cut ",sum(is.Cut.Gene)," genes (", signif(100*sum(is.Cut.Gene)/length(is.Cut.Gene),3), " %) with more than ", prop.failed*100, " % of samples falling below ",count.cutoff," TPM, leaving ", dim(genefilter.eSet)[1]))
  }
  
  return(!is.Cut.Gene)
}
