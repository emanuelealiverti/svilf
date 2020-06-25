margs = commandArgs(trailingOnly=TRUE)
set.seed(3103)
suppressMessages(require(igraph,quietly = T))
dat_gen = function(n_nodes) {
	#NUMBER OF NODES
	#V = 100 
	#V = 500 
	V = n_nodes
	#V = 2500
	dir.create(path = paste0("V",V),showWarnings = F)
	dir.create(path = paste0("V",V,"/RESULTS"),showWarnings = F)

	#+++++++++++++++++++++
	# Latent Factor Model
	#+++++++++++++++++++++
	H = 2 
	Y = matrix(0,V,V)
	Lt = lower.tri(Y)
	low_tr_id = which(Lt,arr.ind =T)
	Z = matrix(rnorm(V*H),V) 
	L = diag(c(3,-3))

	Y[Lt] = rbinom(sum(Lt), 1, plogis(Z %*% L %*% t(Z) -3)[Lt] )

	y_el = get.edgelist(graph.adjacency(Y))
	(fname = paste0("V",V,"/Scen1_V",V,".RData"))
	(fname_el = paste0("V",V,"/Scen1_V",V,"el.RData"))
	save(y_el, file = fname_el)

	Y = Y + t(Y)
	save(Y, file = fname)

	#++++++++++++++++++
	#  Latent Distance
	#++++++++++++++++++

	H = 2 
	Y = matrix(0,V,V)
	Lt = lower.tri(Y)
	low_tr_id = which(Lt,arr.ind =T)
	Z = matrix(rnorm(V*H),V) 

	Y[Lt] = rbinom(sum(Lt), 1, plogis( c(dist(Z)) - 3))

	y_el = get.edgelist(graph.adjacency(Y))
	(fname = paste0("V",V,"/Scen2_V",V,".RData"))
	(fname_el = paste0("V",V,"/Scen2_V",V,"el.RData"))
	save(y_el, file = fname_el)

	Y = Y + t(Y)
	save(Y, file = fname)

	#++++++++++++++++
	# Stochastic BM
	#++++++++++++++++
	H = 2
	id = sample(1:H,size = V,rep=T)
	mm = model.matrix(~-1+factor(id))
	K = matrix(c(.6,.2,.2,.6),2,2)
	Y = matrix(0,V,V)
	Y[Lt] = rbinom(sum(Lt), 1, (mm %*% K %*% t(mm))[Lt] )

	y_el = get.edgelist(graph.adjacency(Y))
	(fname = paste0("V",V,"/Scen3_V",V,".RData"))
	(fname_el = paste0("V",V,"/Scen3_V",V,"el.RData"))
	save(y_el, file = fname_el)

	Y = Y + t(Y)
	save(Y, file = fname)
}
for(v in margs) dat_gen(as.numeric(v))
#for (v in c(100,500,1000, 5000)) dat_gen(v) 
#for (v in length(margs)) dat_gen(margs[v]) 
