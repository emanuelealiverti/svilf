#ifndef SVILF_H 
#define SVILF_H

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <sys/ioctl.h>

/*
 * get terminal width
 */
int consL()
{
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_col; 
}

/*
 * compute plogis and normalise
 */
arma::vec norm_ilogis(arma::vec x, double inter) {
	arma::vec pr = 1 + exp(-inter-x);
	return normalise(pow(pr, -1), 1);
}

/*
 * compute prnorm and normalise
 */
arma::vec norm_iprobit(arma::vec x, double inter) {
	arma::vec pr(x.n_elem);
	for(int i = 0; i < x.n_elem;i++) {
		pr(i) = R::pnorm(inter+x(i), 0.0, 1.0, 1, 0);
	}
	return normalise(pr, 1);
}

arma::uvec arma_sample(int & N, arma::vec & pvec) {
    // count number of elements in pvec
    arma::uword K = pvec.n_elem;
    // create integers 0 to (K-1) to sample from
    arma::uvec opts = arma::linspace<arma::uvec>(0L, K - 1L, K);
    // sample integer
    return arma::conv_to<arma::uvec>::from(Rcpp::RcppArmadillo::sample(opts, N, true, pvec));
}

#endif
