data {
	int<lower=1> V; // nodes
	int<lower=2> H; // latent space
	int<lower=0> Y[V,V]; //data
	int<lower=0> ntrials;
}

parameters{
	matrix[V,H] Z;
	real alpha;
}

model {
	real nu;
	for(v in 1:V){
		Z[v,] ~ multi_normal( rep_vector(0.0,H), diag_matrix(rep_vector(1.0,H)));
	}
	alpha ~ normal(0.0, 2.0);

	// loop over observartion
	for(u in 2:V)
	{
		for(v in 1:u)
		{
			nu = alpha + dot_product(to_vector(Z[u, ]), Z[v,]);
			Y[u,v] ~ binomial(ntrials, inv_logit(nu));
		}
	}
}
