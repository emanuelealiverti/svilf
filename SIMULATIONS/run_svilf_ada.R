#++++++++++++++++++++++++
# read arguments from CLI
#++++++++++++++++++++++++
seed = 27
margs = commandArgs(trailingOnly=TRUE)
nodes = margs[1]
scen = margs[2]
(data_fold = paste0("V",nodes))
(res_fold = paste0(data_fold,"/RESULTS/"))

(file_name = paste0(data_fold,"/Scen",scen,"_",data_fold,"el.RData"))

(res_name = paste0(res_fold,"/Scen",scen,"_svilf_ada.RData"))

require(svilf,lib.loc = "COMPILE/lib")
dir.create(res_fold,showWarnings = F)
load(file_name)
H = 4
set.seed(seed)
tmp = svilf(y_el-1, H = H, intercept=T,get_samples = 2500, prop = 2.0,
	    opts = svilf_options(print_each = 2, tol = 1e-5, 
				 maxit = 5000,svi.maxit_inner =100,sample_adaptive = T))
save(tmp, file=res_name)
#++++++++++++++++++++
# Compute predictions
#++++++++++++++++++++
#niter = 2500 
#V = dim(tmp$E_W)[1]
#Lt = lower.tri(matrix(NA,V,V))
#gc()
#res = matrix(NA, sum(Lt), niter)
#z_samp = tmp$E_W
#for(it in 1:niter) {
	#pres = tmp$a + tcrossprod(tmp$z_samp[,,it])[Lt]
	#res[,it] = plogis(pres)
#}
#pr_hat = rowMeans(res)
#res_name
#save(pr_hat, file = res_name)
