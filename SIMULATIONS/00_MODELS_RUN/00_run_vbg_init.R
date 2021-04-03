# Change directory to work in the top-level one
suppressMessages(suppressWarnings(library(R.utils)))
margs =  commandArgs(asValues=T)
seed  =  margs$seed
nodes =  margs$nodes
scen  =  margs$scen

if(any(sapply(list(seed,nodes,scen), is.null))) stop("USAGE: Rscript -seed X -nodes V -scen S")

(data_fold = sprintf("V%s/seed%s",nodes,seed))

(res_fold = sprintf("%s/RESULTS", data_fold, seed))
(input_name = sprintf("%s/Scen%s_V%s.RData", data_fold, scen, nodes))
(output_name = paste0(res_fold,"/Scen",scen,"_vbg.RData"))


require(rstan)
#dir.create(res_fold,showWarnings = F)
#tmp = stan_model("COMPILE/factor.stan")
if(!file.exists(input_name)) system(sprintf("Rscript ./00_generate_data.R %s %s",seed, nodes))
load('./00_COMPILATION/factor.RData')
load(input_name)

H = 4
Y = as.matrix(Y)
Lt = lower.tri(Y)
ss = svd(Y,nu = H,nv = H)
init = ss$u
rm(ss); gc();
set.seed(seed)

start_at = Sys.time()
fit = vb(tmp, data = list(V = NROW(Y), H = H, Y = Y, ntrials = 1), 
	 iter = 5000, init=list("Z"=init))
pars = extract(fit)
pars$Z = apply(pars$Z, c(2,3), mean) 
pars$alpha = mean(pars$alpha)
pars$lp__ = NULL
end_at = difftime(Sys.time(), start_at, units = 'secs')
mem_all = pryr::mem_used()
pars$time = end_at
pars$memory = mem_all

save(pars, file = output_name)
