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
(output_name = paste0(res_fold,"/Scen",scen,"_lvm4net.RData"))

dir.create(res_fold,showWarnings = F)
load(input_name)


# Fails very often
library(lvm4net)
H = 4

start_at = Sys.time()
set.seed(seed)
tmp = lsm(Y[,], D = H, randomZ = F,nstart = 3,tol = 1e-2)

end_at = difftime(Sys.time(), start_at, units = 'secs')
mem_all = pryr::mem_used()
tmp$time = end_at
tmp$memory = mem_all

save(tmp, file=output_name)
