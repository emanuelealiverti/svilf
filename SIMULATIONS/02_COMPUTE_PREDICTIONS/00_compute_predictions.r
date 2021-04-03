#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compute posterior preditive via monte carlo integration
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
suppressMessages(suppressWarnings(library(R.utils)))
#margs = list(seed=1,nodes=100,scen=1,method='svilf')
margs =  commandArgs(asValues=T)
seed  =  margs$seed
nodes =  margs$nodes
scen  =  margs$scen
method = margs$method

if(any(sapply(list(method, seed,nodes,scen), is.null))) stop("USAGE: Rscript -seed X -method A -nodes V -scen S")

(data_fold = sprintf("V%s/seed%s",nodes,seed))
(input_name = sprintf("%s/Scen%s_V%s.RData", data_fold, scen, nodes))
load(input_name)
y_truth = Y[lower.tri(Y)]
(res_fold = sprintf("%s/RESULTS", data_fold, seed))


(pred_fold = sprintf("OUTPUT/%s/", data_fold, seed))
dir.create(paste0(pred_fold, "/TAB/"),showWarnings = F,recursive = T)
(res_file = paste0(pred_fold, "res_metr", scen,".txt"))
(res_perf = paste0(pred_fold, "perf_metr", scen,".txt"))


(res_name = paste0(res_fold,"/Scen",scen,"_", method, ".RData"))


metr_list = c("auc", "gini","acc","recallTR","precisionTR")

require(ModelMetrics)
acc         =  function (actual, predicted) { sum(diag(table(actual, predicted > mean(actual))))/length(actual)}
precisionTR =  function (actual, predicted) { precision(actual, predicted, cutoff = mean(actual))}
recallTR    =  function (actual, predicted) {recall(actual, predicted, cutoff = mean(actual))}






#++++++++++++++++++++
# Compute and evaluate predictions
#++++++++++++++++++++
cat(getwd(), '\n')
load(res_name)

if(any(c("svilf","svilf_ada") == method)) {
	cat("computing predictions for ", method, '\n')
	V = dim(tmp$E_W)[1]
	Lt = lower.tri(matrix(NA,V,V))
	z_samp = tmp$a + tcrossprod(tmp$E_W)[Lt]
	pr_hat = plogis(z_samp)

	used_time = as.numeric(as.difftime(tmp$time,units="secs"))
	used_memory = as.numeric(tmp$memory) / 2^20

} else if (("vbg" == method)) {
	cat("computing predictions for ", method, '\n')
	V = dim(pars$Z)[1]
	Lt = lower.tri(matrix(NA,V,V))
	z_samp = pars$alpha + tcrossprod(pars$Z)[Lt]
	pr_hat = plogis(z_samp)

	used_time = as.numeric(as.difftime(pars$time,units="secs"))
	used_memory = as.numeric(pars$memory) / 2^20

}  else if ("vblpcm" == method) {
	V = tmp$N
	Lt = lower.tri(matrix(NA,V,V))
	suppressWarnings(suppressMessages(library(VBLPCM)))
	cat("computing predictions for ", method, '\n')
	pr_hat = predict.vblpcm(tmp)[Lt]

	used_time = as.numeric(as.difftime(tmp$time,units="secs"))
	used_memory = as.numeric(tmp$memory) / 2^20

}  else if ("lvm4net" == method) {
	cat("computing predictions for ", method, '\n')
	suppressWarnings(suppressMessages(library(lvm4net)))
	pr_hat =(plogis(tmp$xiT + c((dist(tmp$lsmEZ))^2)))

	used_time = as.numeric(as.difftime(tmp$time,units="secs"))
	used_memory = as.numeric(tmp$memory) / 2^20
}

r = sapply(metr_list, function(x) get(x)(actual = y_truth,predicted = pr_hat))
res = t(as.matrix(r))
rownames(res) = method

perf = t(as.matrix(c('time' = used_time, 'mem' = used_memory)))
rownames(perf) = method

res_file
write.table(res, file = res_file,append = T,col.names = !file.exists(res_file))
write.table(perf, file = res_perf,append = T,col.names = !file.exists(res_perf))
