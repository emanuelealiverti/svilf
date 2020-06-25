#+++++++++++++++++++++++++++++++++++
# Load predictions and compare with truth
#+++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++
# read arguments from CLI
#++++++++++++++++++++++++
margs  =  commandArgs(trailingOnly=TRUE)
nodes  =  margs[1]
scen   =  margs[2]

(data_fold = paste0("V",nodes))
(file_name = paste0(data_fold,"/Scen",scen,"_",data_fold,".RData"))
(res_fold = paste0(data_fold,"/RESULTS/"))
(pred_fold = paste0(data_fold,"/PREDICTIONS/"))
dir.create(paste0(pred_fold, "TAB/"))
(res_tab = paste0(pred_fold, "TAB/res_metr", scen,".txt"))



load(file_name)
y_truth = Y[lower.tri(Y)]
metr_list = c("auc", "gini","acc","recallTR","precisionTR")

require(ModelMetrics)
acc = function (actual, predicted) { sum(diag(table(actual, predicted > mean(actual))))/length(actual)}
precisionTR = function (actual, predicted) { precision(actual, predicted, cutoff = mean(actual))} 
recallTR = function (actual, predicted) {recall(actual, predicted, cutoff = mean(actual))} 


f = function(mm) {
	(res_name = paste0(res_fold,"/Scen",scen,"_",mm, ".RData"))
	(pred_name = paste0(pred_fold,"/Scen",scen,"_",mm, ".RData"))
	load(pred_name)
	r = sapply(metr_list, function(x) get(x)(actual = y_truth,predicted = pr_hat))
	names(r) = metr_list
	return(r)
}

(res = sapply(list("svilf", "svilf_ada", "vbg"), function(l) f(l)))
colnames(res) = c("svilf", "svilf_ada", "vbg") 
write.table(res, file = res_tab)

