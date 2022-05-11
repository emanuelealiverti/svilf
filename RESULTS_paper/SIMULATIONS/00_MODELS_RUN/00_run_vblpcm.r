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
(output_name = paste0(res_fold,"/Scen",scen,"_vblpcm.RData"))


dir.create(res_fold,showWarnings = F)
if(!file.exists(input_name)) system(sprintf("Rscript ./00_generate_data.R %s %s",seed, nodes))
load(input_name)
library(VBLPCM)
H = 4
Y = network(y_el)
set.seed(seed)
start_at = Sys.time()
vbstart = vblpcmstart(Y, G = 1, d = H)
tmp  = vblpcmfit(vbstart, STEPS = 1000)

end_at = difftime(Sys.time(), start_at, units = 'secs')
mem_all = pryr::mem_used()
tmp$time = end_at
tmp$memory = mem_all

save(tmp, file=output_name)
