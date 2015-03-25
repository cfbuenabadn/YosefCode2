library(plotrix)
library(aqfig)

YosefPlot = function(x,gf.vec,gene.set, draw.lines = T,out.dir = NULL, plot.name = "yosefplot.pdf"){
  
  # Create output directory, if necessary
  if (!is.null(out.dir) && !file.exists(out.dir)){
    dir.create(out.dir)
  }
  
  NUM_BINS = 20
  require(plotrix)
  require(aqfig)
  y = x[gf.vec,][gsub("_variant[[:digit:]]*","",rownames(x)[gf.vec]) %in% gene.set,]
  y = log10(y+1)
  y = y[order(rowMeans(y)),]
  
  max_log = max(y)
  
  yosef_h = NULL
  for (i in 1:dim(y)[1]){
    yosef_h = rbind(yosef_h,table(cut(unlist(y[i,]),breaks = (0:(NUM_BINS))*max_log/NUM_BINS,include.lowest = T,right = T)))
  }
  
  pdf(paste0(out.dir,"/",plot.name))
  
  image(t(yosef_h), xlab = "log10(TPM + 1)",ylab = "Mean-Sorted Transcripts", col = color.scale(x = (1:100)/100,extremes = c("cyan","red")), xaxt = 'n', yaxt = 'n')
  axis(side = 1, at = -1/(NUM_BINS-1) + (NUM_BINS*(1:floor(max_log))/max_log)/(NUM_BINS-1), labels = 1:floor(max_log))
  vertical.image.legend(zlim = c(0,100), col = color.scale(x = (1:100)/100,extremes = c("cyan","red")))
  
  if(draw.lines){
    lines(-1/(NUM_BINS-1) + (NUM_BINS*(apply(y,1,mean))/max_log)/(NUM_BINS-1),(1:dim(y)[1])/dim(y)[1], lty = 1)
    lines(-1/(NUM_BINS-1) + (NUM_BINS*(apply(y,1,mean) + apply(y,1,sd) )/max_log)/(NUM_BINS-1),(1:dim(y)[1])/dim(y)[1], lty = 2)
    lines(-1/(NUM_BINS-1) + (NUM_BINS*(apply(y,1,mean) - apply(y,1,sd) )/max_log)/(NUM_BINS-1),(1:dim(y)[1])/dim(y)[1], lty = 2)
  }
  
  dev.off()
  
  
}
