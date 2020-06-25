#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compute posterior preditive via monte carlo integration
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# (this version avoids massively large objects leading to memory issues)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
	res = matrix(NA, sum(Lt), 1)
	pres = tmp$a + tcrossprod(tmp$z_samp[,,1])[Lt]
	res = plogis(pres)
	z_samp = tmp$E_W
	for(it in 2:niter) {
		pres = tmp$a + tcrossprod(tmp$z_samp[,,it])[Lt]
		res = res +  plogis(pres)
		cat(it,"\n")
	}
} else if (("vbg" == method)) {
	V = dim(pars$Z)[2]
	Lt = lower.tri(matrix(NA,V,V))
	gc()
	res = matrix(NA, sum(Lt), 1)
	pres = pars$alpha[1] + tcrossprod(pars$Z[1,,])[Lt]
	res = plogis(pres)

	for(it in 2:niter) {
		pres = pars$alpha[it] + tcrossprod(pars$Z[it,,])[Lt]
		res = res + plogis(pres)
		cat(it,"\n")
	}

}

pr_hat = res/niter
save(pr_hat, file = pred_name)
