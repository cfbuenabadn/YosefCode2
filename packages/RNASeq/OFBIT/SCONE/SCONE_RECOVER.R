out_dir = "/data/yosef/users/mbcole/WTP63/Second/SCONE_out/data"
flist = paste(out_dir,list.files(out_dir),sep = "/")

# Scoring
scores = NULL
snames = NULL
for(fil in flist){
  load(fil)
  nom = scone_out$method
  sc = scone_out$evaluation
  if(!is.na(sc)[1]){
    print(nom)
    snames = c(snames,nom)
    scores = rbind(scores,sc$scores)
  }else{
   # print(nom)
  }
}
rownames(scores) = snames
escores = t(t(scores)*c(-1,-1,1,1,1,-1,1))
colnames(escores) = paste(colnames(scores),c("_INV","")[1.5 + c(-1,-1,1,1,1,-1,1)/2],sep = "")

# Score Heatmaps
require(gplots)
heatmap.2(t(apply(scores ,2,rank)[order(rowSums(apply(escores ,2,rank))),]),density.info = 'n',trace = 'n',key.xlab = "Rank",key.title = NA,
          Colv = NA,col = colorRampPalette(c("purple","black","yellow"))(100),
          margins = c(10,10),cexRow = .75,cexCol = .75)
title(main = "Olfactory p63KO Data: Factor-Free Normalization",cex.lab=0.5)

heatmap.2(t(apply(scores ,2,rank)[order(apply((apply(escores ,2,rank)),1,min)),]),density.info = 'n',trace = 'n',key.xlab = "Rank",key.title = NA,
          Colv = NA,col = colorRampPalette(c("purple","black","yellow"))(100),
          margins = c(10,10),cexRow = .75,cexCol = .75)
title(main = "Olfactory p63KO Data: Factor-Free Normalization",cex.lab=0.5)

heatmap.2(t(apply(scores ,2,rank)[order(apply((apply(escores ,2,rank)),1,min)),]),density.info = 'n',trace = 'n',key.xlab = "Rank",key.title = NA,
          col = colorRampPalette(c("purple","black","yellow"))(100),
          margins = c(10,10),cexRow = .75,cexCol = .75)
title(main = "Olfactory p63KO Data: Factor-Free Normalization",cex.lab=0.5)

heatmap.2(t(apply(scores ,2,rank)[order(-apply((apply(escores ,2,rank)),1,min)),][1:10,]),density.info = 'n',trace = 'n',key.xlab = "Score",key.title = NA,
          Colv = NA,col = colorRampPalette(c("purple","black","yellow"))(100),
          margins = c(15,10),cexRow = .75,cexCol = .5)
title(main = "Top 10 SCONE: Olfactory p63KO Data (8714 Genes)",cex.lab=0.5)


load("/data/yosef/users/mbcole/WTP63/Second/SCONE_out/evaluate_material.Rdata")

load("/data/yosef/users/mbcole/WTP63/Second/SCONE_out/data/NOIMPUTE_FQP_NOWEIGHT_NOBIO_NOBATCH_Q_2.Rdata")
plot(scone_out$evaluation$pc_val, xlab = "wPC1", ylab = "wPC2",main = "Top SCONE: Olfactory p63KO Data (8714 Genes)")
load("/data/yosef/users/mbcole/WTP63/Second/SCONE_test_405/data/NOIMPUTE_FQP_NOWEIGHT_NOBIO_NOBATCH_Q_1.Rdata")
plot(scone_out$evaluation$pc_val, xlab = "wPC1", ylab = "wPC2",main = "Top SCONE: Olfactory p63KO Data (12145 Genes)")

dim(scone_out$exp_mat)
wpc = wPCA(log(scone_out$exp_mat+1),1 - evaluate_material$fnr_pi,nu = 5)
plot(wpc$x[,c(,4)])
