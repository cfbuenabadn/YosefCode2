rm(list=ls())
load("~/archive/users/allonwag/OFBIT_scoreLeafError/normreport2.txthello.Rdata")

tf.vec = TFilter(e=e, type="TPM", method="strong")
tf.vec = tf.vec & (apply(e,1,sd) > SD_EPSILON) # Only consider variable genes for PCA

tf.vec = T
EPSILON = 1
epc = prcomp(t(log(e[tf.vec,]+EPSILON)),center = T,scale. = T)
