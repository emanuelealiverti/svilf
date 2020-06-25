#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compute posterior preditive via monte carlo integration
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
margs  =  commandArgs(trailingOnly=TRUE)
#Dataset
ww =  margs[1]
#method
mm =  margs[2]
require(svilf, lib = "./../lib/")
#++++++++++++++++++++
# Compute predictions
#++++++++++++++++++++
load(paste0("RESULTS/",ww,"/res_",mm,".RData"))
( pred_name = paste0("RESULTS/",ww,"/latent_mean",mm,".RData") )
W = tmp$E_W
save(W, file = pred_name)
