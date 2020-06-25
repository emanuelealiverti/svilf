svilf_options = function(maxit = 100L, print_each = 10L, tol = 1e-5, 
			 svi.tau = 1.0, svi.kappa = 0.75, svi.maxit_inner = 300L, sample_adaptive = FALSE) {
	if(maxit < 1) stop("maxit should be at least 1")
	svi.list = list(tau = svi.tau, kappa = svi.kappa, maxit_inner = svi.maxit_inner)
	return(list(maxit = as.integer(maxit), 
		    sample_adaptive = sample_adaptive,
		    print_each = as.integer(print_each), 
		    tol = tol,
		    svi = svi.list))
}

print_options = function(H, intercept,link, opts){

	set_string = sprintf(paste("Running svilf with H = %g,", link,"link function and", ifelse(intercept,"with","without"), "intercept\n"),H)
	opt_string = sprintf("Options provided: maxit %g, printing each %g iter, tol %g.\n SVI options: tau %g, kappa %g, inner iterations %g. \n",
			     opts$maxit, opts$print_each, opts$tol,
			     opts$svi$tau, opts$svi$kappa, opts$svi$maxit_inner)
	cat(set_string, opt_string)
}

svilf = function(e, H = 2, intercept=T, eigen_init = F, prop = 2.0, link = "logit", get_samples = 0, opts = svilf_options()) {
	if(NCOL(e) != 2) stop("e must be an edgelist (n_edges x 2 matrix or data.frame)")
	if(H <= 1) stop("H must be at least 2")
	if(prop < 1) stop("prop should be at least 1")
	#cat(" Options provided : \n")
	print_options(H=H,intercept=intercept,link=link,opts=opts)

	gr = igraph::graph_from_edgelist(el = as.matrix(e)+1, directed = F)
	null_deg = which(igraph::degree(gr) == 0)
	gr = igraph::delete_vertices(graph = gr,v = null_deg)
	Y = igraph::get.adjacency(gr,  sparse = T)

	#out = svilf_interal(Y, H, prop , maxit = opts$maxit, opts$print_each, opts$tol)
	if(link == "logit") {
		out = svilf_interal_logit(Y, H, prop = prop,sample_adaptive = opts$sample_adaptive, 
					  intercept = intercept, 
					   eigen_init = eigen_init,
					  get_samples = get_samples,
					  tau = opts$svi$tau, kappa = opts$svi$kappa, 
					  maxit = opts$maxit, maxit_inner = opts$svi$maxit_inner, 
					  print_each = opts$print_each, tol = opts$tol)
	} else if (link == "probit") {
		out = svilf_interal_probit(Y, H, prop = prop,sample_adaptive = opts$sample_adaptive, 
					   intercept = intercept, 
					   eigen_init = eigen_init,
					  get_samples = get_samples,
					   tau = opts$svi$tau, kappa = opts$svi$kappa, 
					   maxit = opts$maxit, maxit_inner = opts$svi$maxit_inner, 
					   print_each = opts$print_each, tol = opts$tol)
	} else stop("Link function should be either logit or probit")

	return(out)
}
