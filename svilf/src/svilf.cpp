#include "svilf.h"

using namespace Rcpp;
using namespace arma;


/*
 * SECOND VERSION
 *
 */
// [[Rcpp::export]]
Rcpp::List svilf_interal_logit(const arma::sp_mat Y,
		const int H,
		const double prop = 2.0,
		const bool sample_adaptive = false,
		const bool intercept = true,
		const bool eigen_init= false,
		const int get_samples = 0,
		const double tau = 1.0,
		const double kappa = 0.75,
		const int maxit = 100,
		const int maxit_inner = 10,
		const int print_each = 10,
		const double tol = 1e-5) {
	/*
	 * Get quanties and set structures
	 */
	int V = Y.n_rows;
	arma::mat E_W(V, H, arma::fill::randn);
	arma::mat E_W_old(V, H, arma::fill::randn);


	arma::mat E_W_W(H,H, arma::fill::randn);
	arma::mat curr_cov(H,H, arma::fill::eye);
	arma::mat curr_eta2(H,H, arma::fill::eye);
	arma::vec curr_mean(H, arma::fill::randn);
	arma::vec curr_eta1(H, arma::fill::randn);
	arma::vec old_mean(H);


	arma::cube E_W2(H,H, V);
	arma::cube S_W(H,H, V);
	//std::vector<arma::mat> S_W(V);

	double mm = Y.n_nonzero;
	double alpha = 0.0;
	if(intercept) {
		alpha = log(mm) - log( V * (V-1) - mm);
	}
	double W2 = 0.0;
	arma::vec diff(maxit);
	double curr_diff = 99;
	double global_diff = 99;
	int curr_j;
	double N_n, rho;
	double w = 1.0;

	for(int i = 0; i < V; i++) {
		S_W.slice(i) = arma::eye(H,H);
		E_W_W = E_W.row(i).t() * E_W.row(i);
		E_W_W += S_W.slice(i);
		E_W2.slice(i) = E_W_W;
	}



	/*
	 * for printing
	 */
	int dotsL = consL();
	dotsL = ((50) > (dotsL) ? (dotsL) : (50));
	Rcpp::String star("*");
	Rcpp::String bar("=");
	for(int k = 0; k < dotsL-1; k++) {
		star += "*";
		bar += "=";
	}

	if(eigen_init) {
		Rcpp::Rcout << star.get_cstring() << std::endl;
		Rcpp::Rcout << " Initialising at eigenvalue " << std::endl;
		Rcpp::Rcout << star.get_cstring() << std::endl;
		arma::vec eigval;
		eigs_sym(eigval, E_W, Y, H);
	}
	/*
	 * ****************************************************
	 * MAIN LOOP
	 * ****************************************************
	 */

	Rcpp::Rcout << star.get_cstring() << std::endl;
	Rcpp::Rcout << " V is " << V <<  std::endl;
	Rcpp::Rcout << star.get_cstring() << std::endl;

	int current_iter = 1;
	int inner_iter = 1;
	//construct bar

	clock_t begin = clock();
	while(current_iter < maxit & global_diff > tol) {
		/*
		 * Update in random order
		 */
		E_W_old = E_W;
		arma::vec id = arma::shuffle(arma::regspace(0, V-1));
		for(int j = 0; j < V ; j++){
			int curr_v = id[j];

			/*
			 * Number of nonzero elements (nyi) + number of samped zeros (nyo)
			 * nyo is the minimum between V-nyi-1 (the zeros) and prop*nyi
			 */
			int nyi = accu(Y.row(curr_v));
			int nyo = ((V - 1 - nyi) < std::floor(prop*nyi)) ? (V - 1 - nyi) : std::floor(prop*nyi);
			int totS = nyi + nyo;
			N_n = ( V - 1 - nyi ) / nyo;

			//Rcpp::Rcout << "set zeros" << std::endl;
			/*
			 * Now we create the vector of sampled indexes as follows:
			 * Take all the indexes associated with ones
			 * Then sample nyo elements among the indexes of the zeros
			 */
			arma::uvec id_ones = find(Y.row(curr_v));
			arma::uvec id_onesL = id_ones;
			// trick to avoid to pick curr_v
			id_onesL.resize(nyi+1);
			id_onesL(nyi) = curr_v;

			/*
			 * Reset working quantities
			 *
			 */
			inner_iter = 1;
			curr_diff = 99.0;
			curr_eta2 = S_W.slice(curr_v).i();
			curr_eta1 = curr_eta2 * E_W.row(curr_v).t();
			curr_mean.zeros();

			//
			//while(inner_iter < maxit_inner & curr_diff > tol) {
			arma::uvec id_zeros = arma::regspace<arma::uvec>(0,V-1);
			// Require Armadillo > 0.9800
			id_zeros.shed_rows(id_onesL);
			arma::uvec red_id_zeros(nyo);
			arma::uvec curr_ids(nyo);
			if(sample_adaptive == true) {
				arma::vec pred = E_W * E_W.row(curr_v).t();
				arma::vec samp_pr = norm_ilogis(pred.elem(id_zeros), alpha);
				curr_ids =  arma_sample(nyo, samp_pr);
				red_id_zeros = id_zeros.elem(curr_ids);
				N_n = pow( accu(samp_pr.elem(curr_ids)), -1.0);
			} else {
				//uniform sampling (faster)
			//Rcpp::Rcout << "sampling" << std::endl;
				arma::uvec id_red = arma::regspace<arma::uvec>(0, nyo - 1);
				red_id_zeros = id_zeros.elem(arma::shuffle(id_red));
			}

			arma::uvec idJ = join_cols(id_ones, red_id_zeros);
			// with this we are 100% sure to set the right size
			totS = idJ.n_elem;

			// nested loop. Each E_w must reach the global optimum
			while(inner_iter < maxit_inner & curr_diff > tol) {
			//Rcpp::Rcout << "inner loop" << std::endl;


				arma::vec Z(totS);
				E_W_W.zeros();
				for (int ind = 0; ind < totS;ind++){
			//Rcpp::Rcout << "selecting index" << std::endl;
					curr_j = idJ(ind);
			//Rcpp::Rcout << "index selected" << std::endl;
					W2 =  dot( E_W.row(curr_v), E_W.row(curr_j) ) * 2.0 * alpha;
			//Rcpp::Rcout << "rows picked" << std::endl;
					W2 = W2 * 2.0*alpha;
					W2 += pow(alpha, 2.0);

					W2 += accu( E_W2.slice(curr_v) % E_W2.slice( curr_j ) );
			//Rcpp::Rcout << "slice picked" << std::endl;
					Z(ind) = 0.5 * pow(W2, -0.5) * std::tanh(0.5 * pow(W2, 0.5));
					Z(ind) = std::isfinite(Z(ind)) ? Z(ind) : 0.25;
			//Rcpp::Rcout << "id picked" << std::endl;


					//now we update the covariance function.
					// take care of the fact that when we are accounting for zeros, then
					// we should weight it
					w = (ind < nyi) ? 1.0 : N_n;
					E_W_W += ( w*E_W2.slice(curr_j) * Z(ind) );
				}

			//Rcpp::Rcout << "done inner loop" << std::endl;


				arma::mat E_Wc(H,totS);
				arma::vec Yc(totS);
				for (int ind = 0; ind < totS; ind++){
					w = (ind < nyi) ? 1.0 : N_n;
					E_Wc.col(ind) = w * E_W.row( idJ(ind) ).t();
					Yc(ind) = (ind < nyi) ? 0.5 : -0.5;
				}
			//Rcpp::Rcout << "done weights" << std::endl;

				Yc -= (Z*alpha);

				rho = pow(current_iter*inner_iter + tau, -kappa);

				curr_eta2 = ( 1-rho ) * curr_eta2  + rho * ( E_W_W + arma::eye(H,H) );
				curr_eta1 = ( 1-rho ) * curr_eta1 + rho * (E_Wc * Yc);

				// now we update moments (which are not the etas)
				old_mean = curr_mean;

				curr_cov = curr_eta2.i();
				curr_mean = curr_cov * curr_eta1;
				curr_diff = arma::accu(arma::pow(old_mean - curr_mean, 2)) / (curr_mean.n_elem ) ;


				S_W.slice(curr_v) = curr_cov;
				E_W.row(curr_v) = curr_mean.t();
				E_W2.slice(curr_v) = curr_cov + curr_mean * curr_mean.t();

				Rcpp::checkUserInterrupt();
				inner_iter++;
			}
		}
		global_diff = arma::accu(arma::pow(E_W - E_W_old, 2)) / E_W.n_elem;
		diff(current_iter) = curr_diff;

		if( (current_iter % print_each) == 0) {
			Rcpp::Rcout << bar.get_cstring() << std::endl;
			Rcpp::Rcout << "it " << current_iter << " err " << global_diff << std::endl;
			Rcpp::Rcout << bar.get_cstring() << std::endl;
			if(current_iter == print_each) {
				clock_t end = clock();
				double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
				Rcpp::Rcout << bar.get_cstring() << std::endl;
				Rcpp::Rcout << "CPU time: " << elapsed_secs << " seconds for  " << print_each << " iterations. " << std::endl;
				Rcpp::Rcout << bar.get_cstring() << std::endl;
			}
		}
		current_iter++;
	}
	diff.resize(current_iter);
	Rcpp::List out;
	//out["V"] = V;
	out["S_W"] = S_W;
	out["E_W"] = E_W;
	//out["ew"] = E_W2;
	out["a"] = alpha;
	out["err"] = diff;
	if(get_samples > 1) {
		Rcpp::Rcout << bar.get_cstring() << std::endl;
		Rcpp::Rcout << "Sampling from approx posterior"<< std::endl;
		Rcpp::Rcout << bar.get_cstring() << std::endl;

		// working quantities
		arma::mat curr_z(H,V);
		arma::cube z_samp(V,H,get_samples);
		for(int ss = 0; ss<get_samples; ss++) {
			for(int i = 0; i < V; i ++) {
				curr_z.col(i) = arma::mvnrnd(E_W.row(i).t(), S_W.slice(i));
			}
				z_samp.slice(ss) = curr_z.t();
		}
		out["z_samp"] = z_samp;
	}

	return out;

}


/*
 * PROBIT
 */
// [[Rcpp::export]]
Rcpp::List svilf_interal_probit(const arma::sp_mat Y,
		const int H,
		const double prop = 2.0,
		const bool sample_adaptive = false,
		const bool intercept = true,
		const bool eigen_init= false,
		const int get_samples = 0,
		const double tau = 1.0,
		const double kappa = 0.75,
		const int maxit = 100,
		const int maxit_inner = 10,
		const int print_each = 10,
		const double tol = 1e-5) {
	/*
	 * Get quanties and set structures
	 */
	int V = Y.n_rows;
	arma::mat E_W(V, H, arma::fill::randn);


	arma::mat E_W_old(V, H, arma::fill::randn);
	arma::mat E_W_W(H,H, arma::fill::randn);
	arma::mat curr_cov(H,H, arma::fill::eye);
	arma::mat curr_eta2(H,H, arma::fill::eye);
	arma::vec curr_mean(H, arma::fill::randn);
	arma::vec curr_eta1(H, arma::fill::randn);
	arma::vec old_mean(H);


	arma::cube E_W2(H,H, V);
	arma::cube S_W(H,H, V);
	//std::vector<arma::mat> S_W(V);

	double mm = Y.n_nonzero;
	double Y_tilde;
	double alpha = 0.0;
	if(intercept) {
		alpha = R::qnorm(mm/( V*(V-1) ), 0.0,1.0,1,0);
	}
	double W2 = 0.0;
	arma::vec diff(maxit);
	double curr_diff = 99;
	double global_diff = 99;
	int curr_j;
	double N_n, rho;
	double w = 1.0;

	for(int i = 0; i < V; i++) {
		S_W.slice(i) = arma::eye(H,H);
		E_W_W = E_W.row(i).t() * E_W.row(i);
		E_W_W += S_W.slice(i);
		E_W2.slice(i) = E_W_W;
	}



	/*
	 * for printing
	 */
	int dotsL = consL();
	dotsL = ((50) > (dotsL) ? (dotsL) : (50));
	Rcpp::String star("*");
	Rcpp::String bar("=");
	for(int k = 0; k < dotsL-1; k++) {
		star += "*";
		bar += "=";
	}

	if(eigen_init) {
		Rcpp::Rcout << star.get_cstring() << std::endl;
		Rcpp::Rcout << " Initialising at eigenvalue " << std::endl;
		Rcpp::Rcout << star.get_cstring() << std::endl;
		arma::vec eigval;
		eigs_sym(eigval, E_W, Y, H);
	}

	/*
	 * ****************************************************
	 * MAIN LOOP
	 * ****************************************************
	 */

	Rcpp::Rcout << star.get_cstring() << std::endl;
	Rcpp::Rcout << " V is " << V <<  std::endl;
	Rcpp::Rcout << star.get_cstring() << std::endl;

	int current_iter = 1;
	int inner_iter = 1;
	//construct bar

	clock_t begin = clock();
	while(current_iter < maxit & global_diff > tol) {
		/*
		 * Update in random order
		 */
		E_W_old = E_W;
		arma::vec id = arma::shuffle(arma::regspace(0, V-1));
		for(int j = 0; j < V ; j++){
			int curr_v = id[j];

			/*
			 * Number of nonzero elements (nyi) + number of samped zeros (nyo)
			 * nyo is the minimum between V-nyi-1 (the zeros) and prop*nyi
			 */
			int nyi = accu(Y.row(curr_v));
			int nyo = ((V - 1 - nyi) < std::floor(prop*nyi)) ? (V - 1 - nyi) : std::floor(prop*nyi);
			int totS = nyi + nyo;
			N_n = ( V - 1 - nyi ) / nyo;

			/*
			 * Now we create the vector of sampled indexes as follows:
			 * Take all the indexes associated with ones
			 * Then sample nyo elements among the indexes of the zeros
			 */
			arma::uvec id_ones = find(Y.row(curr_v));
			arma::uvec id_onesL = id_ones;
			// trick to avoid to pick curr_v
			id_onesL.resize(nyi+1);
			id_onesL(nyi) = curr_v;

			/*
			 * Reset working quantities
			 *
			 */
			inner_iter = 1;
			curr_diff = 99.0;
			curr_eta2 = S_W.slice(curr_v).i();
			curr_eta1 = curr_eta2 * E_W.row(curr_v).t();
			curr_mean.zeros();

			//
			//while(inner_iter < maxit_inner & curr_diff > tol) {
			arma::uvec id_zeros = arma::regspace<arma::uvec>(0,V-1);
			// Require Armadillo > 0.9800
			id_zeros.shed_rows(id_onesL);
			arma::uvec red_id_zeros(nyo);
			arma::uvec curr_ids(nyo);
			if(sample_adaptive == true) {
				arma::vec pred = E_W * E_W.row(curr_v).t();
				arma::vec samp_pr = norm_iprobit(pred.elem(id_zeros), alpha);

				curr_ids =  arma_sample(nyo, samp_pr);
				red_id_zeros = id_zeros.elem(curr_ids);
				N_n = pow( accu(samp_pr.elem(curr_ids)), -1.0);

			} else {
				//uniform sampling (faster)
				arma::uvec id_red = arma::regspace<arma::uvec>(0, nyo - 1);
				red_id_zeros = id_zeros.elem(arma::shuffle(id_red));
			}

			arma::uvec idJ = join_cols(id_ones, red_id_zeros);
			// with this we are 100% sure to set the right size
			totS = idJ.n_elem;

			// nested loop. Each E_w must reach the global optimum
			while(inner_iter < maxit_inner & curr_diff > tol) {

				arma::vec Z(totS);
				arma::mat E_Wc(H,totS);
				E_W_W.zeros();
				for (int ind = 0; ind < totS;ind++){
					curr_j = idJ(ind);
					Y_tilde = (Y(curr_v, curr_j) == 0.0) ? -1.0 : 1.0;
					W2 = alpha + dot( E_W.row(curr_v), E_W.row(curr_j) );
					Z(ind) = W2 + Y_tilde * R::dnorm(W2,0.0,1.0,0) / R::pnorm(Y_tilde*W2, 0.0, 1.0, 1, 0);
					Z(ind) = std::isfinite(Z(ind)) ? Z(ind) : 0;

					//now we update the covariance function.
					// take care of the fact that when we are accounting for zeros, then
					// we should weight it
					w = (ind < nyi) ? 1.0 : N_n;
					E_W_W += ( w * E_W2.slice(curr_j) );
					E_Wc.col(ind) = w * E_W.row( curr_j ).t();
				}


				Z -= alpha;

				rho = pow(current_iter*inner_iter + tau, -kappa);

				curr_eta2 = ( 1-rho ) * curr_eta2  + rho * ( E_W_W + arma::eye(H,H) );
				curr_eta1 = ( 1-rho ) * curr_eta1 + rho * (E_Wc * Z);

				// now we update moments (which are not the etas)
				old_mean = curr_mean;

				curr_cov = curr_eta2.i();
				curr_mean = curr_cov * curr_eta1;
				curr_diff = arma::accu(arma::pow(old_mean - curr_mean, 2)) / (curr_mean.n_elem ) ;


				S_W.slice(curr_v) = curr_cov;
				E_W.row(curr_v) = curr_mean.t();
				E_W2.slice(curr_v) = curr_cov + curr_mean * curr_mean.t();

				Rcpp::checkUserInterrupt();
				inner_iter++;
			}
		}
		global_diff = arma::accu(arma::pow(E_W - E_W_old, 2)) / E_W.n_elem;
		diff(current_iter) = curr_diff;

		if( (current_iter % print_each) == 0) {
			Rcpp::Rcout << bar.get_cstring() << std::endl;
			Rcpp::Rcout << "it " << current_iter << " err " << global_diff << std::endl;
			Rcpp::Rcout << bar.get_cstring() << std::endl;
			if(current_iter == print_each) {
				clock_t end = clock();
				double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
				Rcpp::Rcout << bar.get_cstring() << std::endl;
				Rcpp::Rcout << "CPU time: " << elapsed_secs << " seconds for  " << print_each << " iterations. " << std::endl;
				Rcpp::Rcout << bar.get_cstring() << std::endl;
			}
		}
		current_iter++;
	}
	diff.resize(current_iter);
	Rcpp::List out;
	//out["V"] = V;
	out["S_W"] = S_W;
	out["E_W"] = E_W;
	//out["ew"] = E_W2;
	out["a"] = alpha;
	out["err"] = diff;

	if(get_samples > 1) {
		Rcpp::Rcout << bar.get_cstring() << std::endl;
		Rcpp::Rcout << "Sampling from approx posterior"<< std::endl;
		Rcpp::Rcout << bar.get_cstring() << std::endl;

		// working quantities
		arma::mat curr_z(H,V);
		arma::cube z_samp(V,H,get_samples);
		for(int ss = 0; ss<get_samples; ss++) {
			for(int i = 0; i < V; i ++) {
				curr_z.col(i) = arma::mvnrnd(E_W.row(i).t(), S_W.slice(i));
			}
				z_samp.slice(ss) = curr_z.t();
		}
		out["z_samp"] = z_samp;
	}


	return out;

}
/*
   /+++++++++++++++++++++++++++++++++++
   // Utilites for pre/post processing
   /+++++++++++++++++++++++++++++++++++
*/

// [[Rcpp::export]]
arma::vec eff_lowtri(const arma::sp_mat Y) {
	int V = Y.n_rows;
	arma::vec lt(V*(V-1)/2);
	int id = 0;

	for(int u = 0; u < (V-1); u++)
	{
		for(int v = u+1; v < V; v++)
		{
			lt(id) = Y(v,u);
			id++;
		}
	}
	return lt;
}

// [[Rcpp::export]]
arma::vec eff_cross_lt(const arma::mat Y) {
	int V = Y.n_rows;
	arma::vec lt(V*(V-1)/2);
	arma::vec tt(1);
	int id = 0;

	for(int u = 0; u < (V-1); u++)
	{
		for(int v = u+1; v < V; v++)
		{
			tt = Y.row(u) * (Y.row(v).t());
			lt(id) = tt(0);
			id++;
		}
	}
	return lt;
}

// [[Rcpp::export]]
arma::mat id_lt(const arma::vec vv, int V) {
        int n = vv.n_elem;
        arma::mat id(n, 2);
        arma::vec subd = arma::cumsum(arma::regspace(V-1,-1,1));

        int n_row, idc;
        arma::uvec n_col;
        for(int i = 0; i < n; i++)
        {
                n_col = find(subd >= vv(i), 1,"first");
                idc = (n_col(0) == 0) ? 1 : (subd(n_col(0)-1) + 1 );
                n_row = n_col(0) + 1 + vv(i) - idc;
                id(i,0) = n_row + 1;
                id(i,1) = n_col(0)+1;
                if((i % 5000) == 0) {
                        Rcpp::Rcout << i << " out of " << n << std::endl;
                }
        }
        return id;
}

