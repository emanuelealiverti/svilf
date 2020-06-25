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
( pred_name = paste0("RESULTS/",ww,"/pred_mean_",mm,".RData") )
niter = 2500 
V = NROW(tmp$E_W)
pres = res = numeric(V*(V-1)/2)
pres = as.vector(tmp$a + eff_cross_lt(tmp$E_W))
res = plogis(pres)

save(res, file = pred_name)
