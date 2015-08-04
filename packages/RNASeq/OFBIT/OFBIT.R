source("RNASeqNorm.R")
# hk_genes = boolean vector tagging housekeeping genes in rows of e
# prop.q = proportion of variance in preprocessed quality matrix that you wish to retain in producing quality scores (pca)
## ===== OFBIT: OverFitting By Iterative Testing =====
OFBIT = function(e,type,q,techbatch = NULL,biobatch = NULL, hk_genes = NULL, prop.q = 0.80, out.file,plot.dir,
                 binned.norm.methods=c("ComBat","LocScale","GlobScale","AdjEig"),
                 combat.methods=c("yes","no")){
  
  # Make plot directory
  if (file.exists(plot.dir)){
  } else {
    dir.create(plot.dir)
  }
  
  # Produce legend figure: batch and bio
  pdf(paste0(plot.dir,"/legend.pdf"))
  plot(0, type = 'n',bty="n",xaxt = 'n',yaxt = 'n', ylab = "",xlab = "")
  if(!(is.null(biobatch))){
  cols = rainbow(length(levels(as.factor(biobatch))))
  legend("topright",legend = levels(as.factor(biobatch)),col = cols, pch = 16)
  }
  if(!(is.null(techbatch))){
  cols = rainbow(length(levels(as.factor(techbatch))))
  legend("topleft",legend = levels(as.factor(techbatch)),col = cols, pch = 16)
  }
  dev.off()
  
  # Initialize Report
  strout = paste("filt_method","combat_method","scale_method", "quantile_method","qfilt_flag","norm_method","alt_method","bin_method",sep = "|")
  out = cbind(strout,"EPC1_QPC1_pearson","EPC1_QPC1_spearman","tech_batch_KNN_concordance","phenotype_KNN_concordance")
  write.table(out,file = out.file,quote = F,sep = "\t",row.names = F,col.names = F)
  
  # 1) Strong vs weak filtering
  for (filt_method in c("strong","weak")){
    tf.vec = TFilter(e=e, type=type)
    print(filt_method)
    
    # 2) ComBat method: Use or don't use...
    for (combat_method in combat.methods){
      print(combat_method)
      
      if (combat_method == "no" || is.null(techbatch)){
        ce = e
      }else{
        ce = FastComBat(e = e,batch = techbatch,biobatch = biobatch)
      }
      
      # 3) Scale method
      for (scale_method in c(NA,"FQ","UQ","DESeq")){
        
        if(combat_method == "yes" && is.null(techbatch)) break
        
        print(scale_method)
        
        # 4) Quantile method
        if(!is.na(scale_method) && scale_method == "FQ"){
          quantile.methods = c("mean_rank","positive_only")
        }else{
          quantile.methods = NA
        }
        for (quantile_method in quantile.methods){
          print(quantile_method)
          
          if (!is.na(scale_method)){
            sce = EScale(e = ce,tf.vec = tf.vec,method = scale_method,quantile.method = quantile_method)
          }else{
            sce = e
          }
          
          # 5) Filter quality features before computing scores?
          for (qfilt_flag in c(T,F)){
            print(qfilt_flag)
            if (any(is.na(sce))){
              break
            }
            
            # Calculate scores
            score_obj = QPCScores(e = sce,q = q,sig.test = qfilt_flag,tf.vec = tf.vec)
            K = which(cumsum(score_obj$sdev^2/(sum(score_obj$sdev^2))) > prop.q)[1]
            
            #RUV
            if(!is.null(hk_genes)){
              norm_method = "RUVg"
              alt_method = NA
              bin_method = NA
              strout = paste(filt_method,combat_method,scale_method, quantile_method,qfilt_flag,norm_method,alt_method,bin_method,sep = "|")
              print(strout)
              nsce = runRUVg(e = sce,hk_genes,K = K)
              ScoreLeaf(nsce,q,tf.vec = tf.vec,techbatch = techbatch,biobatch = biobatch,leaf.name =  strout,out.file = out.file, plot.dir = plot.dir)
            }
            
            #Regression
            for (norm_method in c("QPCResLoc","QPCResEig","QPCResGlob")){
              alt_method = NA
              bin_method = NA
              strout = paste(filt_method,combat_method,scale_method, quantile_method,qfilt_flag,norm_method,alt_method,bin_method,sep = "|")
              print(strout)
              if(norm_method == "QPCResLoc"){
                nsce = QPCResLoc(e = sce,scores = score_obj$x[,1:K])
              }else if(norm_method == "QPCResEig"){
                nsce = QPCResEig(e = sce,scores = score_obj$x[,1:K])
              }else {
                nsce = QPCResGlob(e = sce,scores = score_obj$x[,1:K])
              }
              ScoreLeaf(nsce,q,tf.vec = tf.vec,techbatch = techbatch,biobatch = biobatch,leaf.name =  strout,out.file = out.file,plot.dir = plot.dir)
            }
            
            #Binned
            for (norm_method in binned.norm.methods){
              
              if (norm_method %in%  c("LocScale","GlobScale")){
                alt.methods = c("median","mean","UQ")
              }else if (norm_method =="AdjEig"){
                alt.methods = c("median","mean")
              }else{
                alt.methods = NA
              }
              
              for (alt_method in alt.methods ){
                
                for (bin_method in c("quantiles","mix_norm")){
                  nsce = sce
                  strout = paste(filt_method,combat_method,scale_method, quantile_method,qfilt_flag,norm_method,alt_method,bin_method,sep = "|")
                  print(strout)
                  
                  for(k in K:1){
                    bins = QBins(score = score_obj$x[,k],method = bin_method)  
                    print(table(bins))                   
                    if(norm_method == "ComBat"){
                      if(any(table(bins) < 2)){
                        print("Disqualified ComBat due to single-sample bin...")
                      }else{
                        nsce = FastComBat(e = nsce,batch = bins,biobatch = biobatch)
                      }
                    }else if (norm_method == "LocScale"){
                      nsce = LocScale(e = nsce,batch = bins,method = alt_method)
                    }else if (norm_method == "GlobScale"){
                      nsce = GlobScale(e = nsce,batch = bins,method = alt_method)
                    }else {
                      nsce = AdjEig(e = nsce,batch = bins,method = alt_method)
                    }
                  }
                  if(norm_method != "ComBat" || any(table(bins) > 1) ){
                    ScoreLeaf(nsce,q,tf.vec = tf.vec,techbatch = techbatch,biobatch = biobatch,leaf.name =  strout,out.file = out.file,plot.dir = plot.dir)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

ScoreLeaf = function(e,q,tf.vec = T,techbatch = NULL,biobatch = NULL, knn = 10, leaf.name, out.file,plot.dir,EPSILON = 1){
  
  epc = prcomp(t(log(e[tf.vec,][apply(e[tf.vec,],1,sd)>0,]+EPSILON)),center = T,scale. = T)
  ppq = PPQual(q)
  qpc = prcomp(ppq)
  cor_pearson = cor(epc$x[,1],qpc$x[,1], method = "pearson")
  cor_spearman = cor(epc$x[,1],qpc$x[,1], method = "spearman")
  
  if(!(is.null(techbatch))){
    
    pdf(paste0(plot.dir,"/",leaf.name,"_epc1vqpc1_techcol.pdf"))
    cols = rainbow(length(levels(as.factor(techbatch))))[as.factor(techbatch)]
    plot(epc$x[,1],qpc$x[,1],pch = 16, col = cols,xlab = "EPC1",ylab = "QPC1")
    dev.off()
    
    pdf(paste0(plot.dir,"/",leaf.name,"_epc1vepc2_techcol.pdf"))
    cols = rainbow(length(levels(as.factor(techbatch))))[as.factor(techbatch)]
    plot(epc$x[,1],epc$x[,2],pch = 16, col = cols,xlab = "EPC1",ylab = "EPC2")
    dev.off()
    
    d = as.matrix(dist(epc$x[,1:2],method = "euclidean"))
    m = NULL
    for(s in 1:dim(d)[1]){
      ds = d[,s]
      ra = rank(ds)
      m[s] = mean(techbatch[(ra <= knn+1) & (ra != 1)] == techbatch[s])
    }
    tech_concord = mean(m)
  }else{
    tech_concord = NA
  }
  
  if(!(is.null(biobatch))){
    pdf(paste0(plot.dir,"/",leaf.name,"_epc1vqpc1_biocol.pdf"))
    cols = rainbow(length(levels(as.factor(biobatch))))[as.factor(biobatch)]
    plot(epc$x[,1],qpc$x[,1],pch = 16, col = cols, xlab = "EPC1",ylab = "QPC1")
    dev.off()
    
    pdf(paste0(plot.dir,"/",leaf.name,"_epc1vepc2_biocol.pdf"))
    cols = rainbow(length(levels(as.factor(biobatch))))[as.factor(biobatch)]
    plot(epc$x[,1],epc$x[,2],pch = 16, col = cols, xlab = "EPC1",ylab = "EPC2")
    dev.off()
    
    d = as.matrix(dist(epc$x[,1:2],method = "euclidean"))
    m = NULL
    for(s in 1:dim(d)[1]){
      ds = d[,s]
      ra = rank(ds)
      m[s] = mean(biobatch[(ra <= knn+1) & (ra != 1)] == biobatch[s])
    }
    bio_concord = mean(m)
  }else{
    bio_concord = NA
  }
  
  
  out = cbind(leaf.name,cor_pearson,cor_spearman,tech_concord,bio_concord)
  write.table(out,file = out.file,append = T,quote = F,sep = "\t",row.names = F,col.names = F)
}