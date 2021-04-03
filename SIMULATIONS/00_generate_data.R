# Generates artificial data over Nseed scenarios w seeds given by 1, 2, ..., Nseeds
library(igraph)
margs = commandArgs(trailingOnly = T)
if(length(margs) < 2) stop("USAGE: Rscript FILE Nseeds V1 V2 ...")
seed = as.numeric(margs[1])

dat_gen = function(seed, n_nodes) {
	#NUMBER OF NODES
	#V = 100 
	#V = 500 
	V = n_nodes
	#V = 2500
	#dir.create(path = paste0("V",V),showWarnings = F,)
	dir.create(path = sprintf("V%s/seed%s/RESULTS", n_nodes, seed),recursive = T, showWarnings = F)

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

	fname = sprintf("V%s/seed%s/Scen1_V%s.RData", n_nodes, seed, n_nodes)
	fname_el = sprintf("V%s/seed%s/Scen1_V%s_el.RData", n_nodes, seed, n_nodes)

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
	fname = sprintf("V%s/seed%s/Scen2_V%s.RData", n_nodes, seed, n_nodes)
	fname_el = sprintf("V%s/seed%s/Scen2_V%s_el.RData", n_nodes, seed, n_nodes)
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

	fname = sprintf("V%s/seed%s/Scen3_V%s.RData", n_nodes, seed, n_nodes)
	fname_el = sprintf("V%s/seed%s/Scen3_V%s_el.RData", n_nodes, seed, n_nodes)

	save(y_el, file = fname_el)

	Y = Y + t(Y)
	save(Y, file = fname)
}
for(v in margs[-1]) dat_gen(seed, as.numeric(v))
