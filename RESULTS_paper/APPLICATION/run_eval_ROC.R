#+++++++++++++++++++++++++++++++++++
# Load predictions and compare with truth
#+++++++++++++++++++++++++++++++++++
margs  =  commandArgs(trailingOnly=TRUE)
dataset =  margs[1]
mm =  margs[2]

require(svilf,lib.loc = "../lib")
require(igraph)
require(Matrix)

filen = dir(paste0("./DATA/",dataset),full.names = T)
dir.create(paste0("RESULTS/",dataset),recursive = T)
( pred_name = paste0("RESULTS/",dataset,"/pred_mean_",mm,".RData") )
(res_roc = paste0("RESULTS/",dataset,"/ROC_",mm, ".RData") )
load(pred_name)
y_el = read.csv(filen)


nn = graph_from_edgelist(el = as.matrix(y_el)+1,directed = F)
mm = as_adjacency_matrix(nn,sparse = T)

y_truth = as.vector(eff_lowtri(mm))

simple_roc <- function(true, scores){
  true <- true[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(true)/sum(true), FPR=cumsum(!true)/sum(!true))
}

r = simple_roc(y_truth,res)
tmp=spline(x=r$TPR,y=r$FPR)
save(tmp, file = res_roc)

