# Change directory to work in the top-level one
suppressMessages(suppressWarnings(library(R.utils)))
margs =  commandArgs(asValues=T)
seed  =  margs$seed
nodes =  margs$nodes
scen  =  margs$scen

if(any(sapply(list(seed,nodes,scen), is.null))) stop("USAGE: Rscript -seed X -nodes V -scen S")

(data_fold = sprintf("V%s/seed%s",nodes,seed))

(res_fold = sprintf("%s/RESULTS", data_fold, seed))
(input_name = sprintf("%s/Scen%s_V%s_el.RData", data_fold, scen, nodes))
(output_name = paste0(res_fold,"/Scen",scen,"_svilf_ada.RData"))

require(svilf,lib.loc = "../lib")
dir.create(res_fold,showWarnings = F)
if(!file.exists(input_name)) system(sprintf("Rscript ./00_generate_data.R %s %s",seed, nodes))
load(input_name)
H = 4
set.seed(seed)

start_at = Sys.time()

tmp = svilf(y_el-1, H = H, intercept=T,get_samples = 1, prop = 2.0,
	    opts = svilf_options(print_each = 2, tol = 1e-5, 
				 maxit = 5000,svi.maxit_inner = 3,sample_adaptive = T))
end_at = difftime(Sys.time(), start_at, units = 'secs')
mem_all = pryr::mem_used()
tmp$time = end_at
tmp$memory = mem_all
save(tmp, file=output_name)
