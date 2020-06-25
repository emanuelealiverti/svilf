#++++++++++++++++++++++++
# read arguments from CLI
#++++++++++++++++++++++++
margs  =  commandArgs(trailingOnly=TRUE)
dataset =  margs[1]
is_ada =  (margs[2] == "svilf_ada")

filen = dir(paste0("./DATA/",dataset),full.names = T)
dir.create(paste0("RESULTS/",dataset),recursive = T)

seed=124
y_el = read.csv(filen)
require(svilf,lib.loc = "../lib")
cat("Running on ", dataset, "dataset", ifelse(is_ada,"with", "without"), "adaptive sampling \n")


H = 5
set.seed(seed)
tmp = svilf(y_el, H = H, intercept=T, eigen_init = T, get_samples = 1, prop = 3.0,
	    opts = svilf_options(print_each = 1, tol = 1e-6, 
				 maxit = 200,svi.maxit_inner = 5,sample_adaptive = is_ada))
(res_name = paste0(paste0("RESULTS/",dataset),"/res",ifelse(is_ada,"_svilf_ada.RData","_svilf.RData")))
save(tmp, file=res_name)
