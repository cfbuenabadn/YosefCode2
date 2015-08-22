source("RNASeqNorm.R")

## ===== OFBIT: OverFitting By Iterative Testing =====
# e = expression matrix (rows = transcripts, cols = samples)
# type = data type (count or TPM)
# techbatch = categorical batch variable - known
# biobatch = categorical biological variable - known - NA is allowed, but NA samples will be disregarded in plots and KNN calculations
# hk_genes = boolean vector tagging housekeeping genes in rows of e
# prop.q = proportion of variance in preprocessed quality matrix that you wish to retain in producing quality scores (qPCA)
# out.file = path to output file (scores etc.)
# tf.vec = transcript filter vector - if not null, will override filtering step as option "user"
# plot.dir = path to plot directory
# ...

OFBIT = function(e,type,q,techbatch = NULL,biobatch = NULL, hk_genes = NULL, prop.q = 0.80, out.file,plot.dir,
                 binned.norm.methods=c("ComBat","LocScale"),
                 combat.methods=c("yes","no"),
                 filtering.methods=c("strong","weak"),
                 scaling.methods = c(NA,"FQ","UQ","DESeq"),
                 tf.vec = NULL,
                 regression.norm.methods=c("ResLoc"),
                 error.file="/dev/null", log.file="/dev/null"){
  
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
  strout = paste("filt_method","scale_method", "zero_method","combat_method","qfilt_flag","norm_method","alt_method","bin_method",sep = "|")
  out = cbind(strout,"EPC12_QPC1_pearson","EPC12_QPC1_spearman","tech_batch_KNN_concordance","phenotype_KNN_concordance","WPC12_QPC1_pearson","WPC12_QPC1_spearman","weighted_tech_batch_KNN_concordance","weighted_phenotype_KNN_concordance")
  write.table(out,file = out.file,quote = F,sep = "\t",row.names = F,col.names = F)
  
  # 1) Strong vs weak transcript filtering
  if(!is.null(tf.vec)){
    filtering.methods = "user"
    user.tf.vec = tf.vec
  }
  for (filt_method in filtering.methods){
    if(filt_method != "user"){
      tf.vec = TFilter(e=e, type=type, method=filt_method)
    }else{
      tf.vec = user.tf.vec
    }
    print(filt_method)
    print(sum(tf.vec))
    
    # 2) Scale method: Scaling data globally / quantile-based
    
    if(sum(rowSums(e[tf.vec,] == 0) == 0) < 2){
      scaling.methods = scaling.methods[scaling.methods != "DESeq"]
      print("< 2 genes are detected in all samples! DESeq will not be run.")
    }
    
    for (scale_method in scaling.methods){
      print(scale_method)
      # 3) Zero-Handling Method
      
        if(!is.na(scale_method) && scale_method == "FQ"){
          zero.methods = c("all","positive")
        }else if(!is.na(scale_method) && scale_method == "UQ"){
          if ( !any(apply(e[tf.vec,],2,quantile,na.rm = T)[4,] == 0) ){
            zero.methods = c("all","positive")
          }else{
            print("UQ = 0 for subset of samples! Only positive method will be used for UQ.")
            zero.methods = c("positive")
          }
        }else{
          zero.methods = NA
        }
        for (zero_method in zero.methods){
          print(zero_method)
          if (!is.na(scale_method)){
            se = EScale(e = e,tf.vec = tf.vec,method = scale_method,zero.method = zero_method )
          }else{
            se = e
          }
          
          # 4) ComBat method: Remove known batch effects prior to other normalization methods?
          if(is.null(techbatch)){
            combat.methods = combat.methods[combat.methods == "no"]
          }
          for (combat_method in combat.methods){
            print(combat_method)
            if (combat_method == "no"){
              ce = se
            }else{
              ce = FastComBat(e = se,batch = techbatch,biobatch = biobatch)
              ce[se == 0] = 0
            }
            
            # 5) Filter quality features before computing scores?
            for (qfilt_flag in c(T,F)){
              print(qfilt_flag)
            
              # Nothing More - Only Once
              if(qfilt_flag){
                norm_method = NA
                alt_method = NA
                bin_method = NA
                strout = paste(filt_method,scale_method, zero_method,combat_method,qfilt_flag,norm_method,alt_method,bin_method,sep = "|")
                print(strout)
                ScoreLeaf(ce,q,tf.vec = tf.vec,techbatch = techbatch,biobatch = biobatch,leaf.name =  strout,out.file = out.file, plot.dir = plot.dir)
              }
              
              # Calculate scores
              score_obj = QPCScores(e = ce,q = q,sig.test = qfilt_flag,tf.vec = tf.vec)            
              K = which(cumsum(score_obj$sdev^2/(sum(score_obj$sdev^2))) > prop.q)[1]
              
              # RUV Regression - Only Once
              if(!is.null(hk_genes) && qfilt_flag){
                
                alt.methods = c("all","positive")
                num_pass_filt = sum(rowSums((ce[hk_genes & tf.vec,] > 0) %*% t(ce[hk_genes & tf.vec ,] > 0) == 0) == 0)
                if (num_pass_filt < 3){
                  print("Not enough weightless control gene comparisons. Positive option invalid.")
                  alt.methods = c("all")
                  max_k = min(5, K)
                }else{
                  max_k = min(5, K, num_pass_filt)
                }
                
                bin_method = NA
              
                for (curK in 1:max_k) {
                  norm_method = paste0("RUVg", "(k=", curK, ")")
                  for(alt_method in alt.methods){
                    strout = paste(filt_method,scale_method, zero_method,combat_method,qfilt_flag,norm_method,alt_method,bin_method,sep = "|")
                    print(strout)
                    ne = YL_RUVg(e = ce,hk_genes & tf.vec,K = curK,zero.method = alt_method)
                    ne[ce == 0] = 0
                    ScoreLeaf(ne,q,tf.vec = tf.vec,techbatch = techbatch,biobatch = biobatch,leaf.name =  strout,out.file = out.file, plot.dir = plot.dir)  
                  }
                }
              }
            
              # Score Regression
              max_k = min(5, K)
              for (norm_method in regression.norm.methods){
                alt.methods = c("all","positive")
                bin_method = NA
                if(norm_method == "ResLoc"){
                  for(alt_method in alt.methods){
                    for (curK in 1:max_k) {
                      norm_method2 = paste0("ResLoc", "(k=", curK, ")")
                      strout = paste(filt_method,scale_method, zero_method,combat_method,qfilt_flag,norm_method2,alt_method,bin_method,sep = "|")
                      print(strout)
                      ne = ResLoc(e = ce,scores = score_obj$x[,1:curK],zero.method = alt_method)
                      ne[ce == 0] = 0
                      ScoreLeaf(ne,q,tf.vec = tf.vec,techbatch = techbatch,biobatch = biobatch,leaf.name =  strout,out.file = out.file,plot.dir = plot.dir)
                    }
                  }
                }
              }
            
            
              # Binned Adjustments for Scores
              for (norm_method in binned.norm.methods){
                if (norm_method %in%  c("LocScale")){
                  alt.methods = c("median_all","mean_all","UQ_all","median_positive","mean_positive","UQ_positive")
                }else{
                  alt.methods = NA
                }
              
              for (alt_method in alt.methods ){
                
                for (bin_method in c("quantiles","mix_norm")){
                  ne = ce
                
                  for(k in 1:K){
                    bins = QBins(score = score_obj$x[,k],method = bin_method)  
                    if(norm_method == "ComBat"){
                      if(any(table(bins) < 2)){
                        print("Disqualified ComBat instance due to single-sample bin...")
                      }else if (length(unique(bins)) == 1){
                        print("Disqualified ComBat instance due to single bin...")
                      }else{
                        ne = FastComBat(e = ne,batch = bins,biobatch = biobatch)
                        ne[ce == 0] = 0
                        norm_method2 = paste0("ComBat", "(k=", k, ")")
                        strout = paste(filt_method,scale_method, zero_method,combat_method,qfilt_flag,norm_method2,alt_method,bin_method,sep = "|")
                        print(strout)
                        ScoreLeaf(ne,q,tf.vec = tf.vec,techbatch = techbatch,biobatch = biobatch,leaf.name =  strout,out.file = out.file,plot.dir = plot.dir)
                      }
                    }else{
                      if (length(unique(bins)) == 1){
                        print("Disqualified LocScale instance due to single bin...")
                      }else{
                        upa = strsplit(alt_method,split = "_")[[1]]
                        ne = LocScale(e = ne,batch = bins,method = upa[1],zero.method = upa[2])
                        ne[ce == 0] = 0
                        norm_method2 = paste0("LocScale", "(k=", k, ")")
                        strout = paste(filt_method,scale_method, zero_method,combat_method,qfilt_flag,norm_method2,alt_method,bin_method,sep = "|")
                        print(strout)
                        ScoreLeaf(ne,q,tf.vec = tf.vec,techbatch = techbatch,biobatch = biobatch,leaf.name =  strout,out.file = out.file,plot.dir = plot.dir)
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
    }
  }


ReportError = function(leaf.name, error.file, error.msg){
  
  out = cbind(leaf.name, error.msg)
  write.table(out,file = error.file,append = T,quote = F,sep = "\t",row.names = F,col.names = F)
  
}

ScoreLeaf = function(e,q,tf.vec = T,techbatch = NULL,biobatch = NULL, knn = 10, leaf.name, out.file,plot.dir,EPSILON = 1){
  
    tf.vec = tf.vec & (apply(e,1,sd) > SD_EPSILON) # Only consider variable genes for PCA
  
    epc = prcomp(t(log(e[tf.vec,]+EPSILON)),center = T,scale. = T)
    wpc = wPCA(log(e[tf.vec,]+EPSILON),w = 0 + (e[tf.vec,] > 0), nu = 2,filt = T,scale. = T)
    
    ppq = PPQual(q)
    qpc = prcomp(ppq)
    cor_pearson = max(abs(cor(epc$x[,1:2],qpc$x[,1], method = "pearson")))
    cor_spearman = max(abs(cor(epc$x[,1:2],qpc$x[,1], method = "spearman")))
    wcor_pearson = max(abs(cor(wpc$x[,1:2],qpc$x[,1], method = "pearson")))
    wcor_spearman = max(abs(cor(wpc$x[,1:2],qpc$x[,1], method = "spearman")))
    
    
    if(!(is.null(techbatch))){
      
      cols = rainbow(length(levels(as.factor(techbatch))))[as.factor(techbatch)]
      
      pdf(paste0(plot.dir,"/",leaf.name,"_epc1vqpc1_techcol.pdf"))
      plot(epc$x[,1],qpc$x[,1],pch = 16, col = cols,xlab = "EPC1",ylab = "QPC1")
      dev.off()
      
      pdf(paste0(plot.dir,"/",leaf.name,"_epc1vepc2_techcol.pdf"))
      plot(epc$x[,1],epc$x[,2],pch = 16, col = cols,xlab = "EPC1",ylab = "EPC2")
      dev.off()
      
      pdf(paste0(plot.dir,"/",leaf.name,"_wpc1vqpc1_techcol.pdf"))
      plot(wpc$x[,1],qpc$x[,1],pch = 16, col = cols,xlab = "EPC1",ylab = "QPC1")
      dev.off()
      
      pdf(paste0(plot.dir,"/",leaf.name,"_wpc1vwpc2_techcol.pdf"))
      plot(wpc$x[,1],wpc$x[,2],pch = 16, col = cols,xlab = "EPC1",ylab = "EPC2")
      dev.off()
      
      d = as.matrix(dist(epc$x[,1:2],method = "euclidean"))
      m = NULL
      for(s in 1:dim(d)[1]){
        ds = d[,s]
        ra = rank(ds)
        m[s] = mean(techbatch[(ra <= knn+1) & (ra != 1)] == techbatch[s])
      }
      tech_concord = mean(m)
      
      d = as.matrix(dist(wpc$x[,1:2],method = "euclidean"))
      m = NULL
      for(s in 1:dim(d)[1]){
        ds = d[,s]
        ra = rank(ds)
        m[s] = mean(techbatch[(ra <= knn+1) & (ra != 1)] == techbatch[s])
      }
      wtech_concord = mean(m)
    }else{
      tech_concord = NA
      wtech_concord = NA
    }
    
    if(!(is.null(biobatch))){
      
      cols = rainbow(length(levels(as.factor(biobatch))))[as.factor(biobatch)]
      cols[is.na(cols)] = "black"
      
      pdf(paste0(plot.dir,"/",leaf.name,"_epc1vqpc1_biocol.pdf"))
      plot(epc$x[,1],qpc$x[,1],pch = 16, col = cols,xlab = "EPC1",ylab = "QPC1")
      dev.off()
      
      pdf(paste0(plot.dir,"/",leaf.name,"_epc1vepc2_biocol.pdf"))
      plot(epc$x[,1],epc$x[,2],pch = 16, col = cols,xlab = "EPC1",ylab = "EPC2")
      dev.off()
      
      pdf(paste0(plot.dir,"/",leaf.name,"_wpc1vqpc1_biocol.pdf"))
      plot(wpc$x[,1],qpc$x[,1],pch = 16, col = cols,xlab = "EPC1",ylab = "QPC1")
      dev.off()
      
      pdf(paste0(plot.dir,"/",leaf.name,"_wpc1vwpc2_biocol.pdf"))
      plot(wpc$x[,1],wpc$x[,2],pch = 16, col = cols,xlab = "EPC1",ylab = "EPC2")
      dev.off()
      
      is.un = is.na(biobatch)
      
      d = as.matrix(dist(epc$x[,1:2][!is.un,],method = "euclidean"))
      m = NULL
      for(s in 1:dim(d)[1]){
        ds = d[,s]
        ra = rank(ds)
        m[s] = mean(biobatch[!is.un][(ra <= knn+1) & (ra != 1)] == biobatch[!is.un][s])
      }
      bio_concord = mean(m)
      
      d = as.matrix(dist(wpc$x[,1:2][!is.un,],method = "euclidean"))
      m = NULL
      for(s in 1:dim(d)[1]){
        ds = d[,s]
        ra = rank(ds)
        m[s] = mean(biobatch[!is.un][(ra <= knn+1) & (ra != 1)] == biobatch[!is.un][s])
      }
      wbio_concord = mean(m)
      
    }else{
      bio_concord = NA
      wbio_concord = NA
    }
    
    out = cbind(leaf.name,cor_pearson,cor_spearman,tech_concord,bio_concord,wcor_pearson,wcor_spearman,wtech_concord,wbio_concord) 
    write.table(out,file = out.file,append = T,quote = F,sep = "\t",row.names = F,col.names = F) 
}