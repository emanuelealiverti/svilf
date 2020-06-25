#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compute posterior preditive via monte carlo integration
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++
# read arguments from CLI
#++++++++++++++++++++++++
margs  =  commandArgs(trailingOnly=TRUE)
nodes  =  margs[1]
scen   =  margs[2]
method =  margs[3]

(data_fold = paste0("V",nodes))
(res_fold = paste0(data_fold,"/RESULTS/"))
(pred_fold = paste0(data_fold,"/PREDICTIONS/"))

dir.create(pred_fold,showWarnings = F)

(res_name = paste0(res_fold,"/Scen",scen,"_", method, ".RData"))
(pred_name = paste0(pred_fold,"/Scen",scen,"_", method, ".RData"))

#++++++++++++++++++++
# Compute predictions
#++++++++++++++++++++
load(res_name)
niter = 2500 

if(any(c("svilf","svilf_ada") == method)) {
	V = dim(tmp$E_W)[1]
	Lt = lower.tri(matrix(NA,V,V))
	gc()
	res = matrix(NA, sum(Lt), niter)
	z_samp = tmp$E_W
	for(it in 1:niter) {
		pres = tmp$a + tcrossprod(tmp$z_samp[,,it])[Lt]
		res[,it] = plogis(pres)
	}
} else if (("vbg" == method)) {
	V = dim(pars$Z)[2]
	Lt = lower.tri(matrix(NA,V,V))
	gc()
	res = matrix(NA, sum(Lt), niter)
	for(it in 1:niter) {
		pres = pars$alpha[it] + tcrossprod(pars$Z[it,,])[Lt]
		res[,it] = plogis(pres)
	}

}

pr_hat = rowMeans(res)
save(pr_hat, file = pred_name)
