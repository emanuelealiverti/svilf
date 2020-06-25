#++++++++++++++++++++++++
# read arguments from CLI
#++++++++++++++++++++++++
seed = 27
margs = commandArgs(trailingOnly=TRUE)
nodes = margs[1]
scen = margs[2]
(data_fold = paste0("V",nodes))
(res_fold = paste0(data_fold,"/RESULTS/"))

(file_name = paste0(data_fold,"/Scen",scen,"_",data_fold,".RData"))

(res_name = paste0(res_fold,"/Scen",scen,"_vbg.RData"))

require(rstan)
dir.create(res_fold,showWarnings = F)
#tmp = stan_model("COMPILE/factor.stan")
load('./COMPILE/factor.RData')
load(file_name)

H = 4
Y = as.matrix(Y)
Lt = lower.tri(Y)
ss = svd(Y,nu = H,nv = H)
init = ss$u
rm(ss); gc();
set.seed(seed)
fit = vb(tmp, data = list(V = NROW(Y), H = H, Y = Y, ntrials = 1), 
	 iter = 5000, 
	 output_samples=2500, init=list("Z"=init))
pars = extract(fit)
save(pars, file = res_name)

#++++++++++++++++++++
# Compute predictions
#++++++++++++++++++++
#niter = length(pars$alpha)
#res = matrix(NA, sum(Lt), niter)
#for(it in 1:niter) {
	#pres = pars$alpha[it] + tcrossprod(pars$Z[it,,])[Lt]
	#res[,it] = plogis(pres)
#}
#pr_hat = rowMeans(res)

#save(pr_hat, file = res_name)
