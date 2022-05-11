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
load(pred_name)
y_el = read.csv(filen)

(res_tab = paste0("RESULTS/",dataset,"/tab_",mm, ".RData") )

nn = graph_from_edgelist(el = as.matrix(y_el)+1,directed = F)
mm = as_adjacency_matrix(nn,sparse = T)

y_truth = as.vector(eff_lowtri(mm))
y_truth = 1*(y_truth != 0)
summary(y_truth)


metr_list = c("auc", "gini","acc","recallTR","precisionTR")

require(ModelMetrics)
acc = function (actual, predicted) { sum(diag(table(actual, predicted > mean(actual))))/length(actual)}
precisionTR = function (actual, predicted) { precision(actual, predicted, cutoff = mean(actual))} 
recallTR = function (actual, predicted) {recall(actual, predicted, cutoff = mean(actual))} 


r = sapply(metr_list, function(x) get(x)(actual = y_truth,predicted = res))
names(r) = metr_list
r
save(r, file = res_tab)

