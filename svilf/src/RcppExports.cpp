// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// svilf_interal_logit
Rcpp::List svilf_interal_logit(const arma::sp_mat Y, const int H, const double prop, const bool sample_adaptive, const bool intercept, const bool eigen_init, const int get_samples, const double tau, const double kappa, const int maxit, const int maxit_inner, const int print_each, const double tol);
RcppExport SEXP _svilf_svilf_interal_logit(SEXP YSEXP, SEXP HSEXP, SEXP propSEXP, SEXP sample_adaptiveSEXP, SEXP interceptSEXP, SEXP eigen_initSEXP, SEXP get_samplesSEXP, SEXP tauSEXP, SEXP kappaSEXP, SEXP maxitSEXP, SEXP maxit_innerSEXP, SEXP print_eachSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type H(HSEXP);
    Rcpp::traits::input_parameter< const double >::type prop(propSEXP);
    Rcpp::traits::input_parameter< const bool >::type sample_adaptive(sample_adaptiveSEXP);
    Rcpp::traits::input_parameter< const bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type eigen_init(eigen_initSEXP);
    Rcpp::traits::input_parameter< const int >::type get_samples(get_samplesSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit_inner(maxit_innerSEXP);
    Rcpp::traits::input_parameter< const int >::type print_each(print_eachSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(svilf_interal_logit(Y, H, prop, sample_adaptive, intercept, eigen_init, get_samples, tau, kappa, maxit, maxit_inner, print_each, tol));
    return rcpp_result_gen;
END_RCPP
}
// svilf_interal_probit
Rcpp::List svilf_interal_probit(const arma::sp_mat Y, const int H, const double prop, const bool sample_adaptive, const bool intercept, const bool eigen_init, const int get_samples, const double tau, const double kappa, const int maxit, const int maxit_inner, const int print_each, const double tol);
RcppExport SEXP _svilf_svilf_interal_probit(SEXP YSEXP, SEXP HSEXP, SEXP propSEXP, SEXP sample_adaptiveSEXP, SEXP interceptSEXP, SEXP eigen_initSEXP, SEXP get_samplesSEXP, SEXP tauSEXP, SEXP kappaSEXP, SEXP maxitSEXP, SEXP maxit_innerSEXP, SEXP print_eachSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const int >::type H(HSEXP);
    Rcpp::traits::input_parameter< const double >::type prop(propSEXP);
    Rcpp::traits::input_parameter< const bool >::type sample_adaptive(sample_adaptiveSEXP);
    Rcpp::traits::input_parameter< const bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const bool >::type eigen_init(eigen_initSEXP);
    Rcpp::traits::input_parameter< const int >::type get_samples(get_samplesSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit_inner(maxit_innerSEXP);
    Rcpp::traits::input_parameter< const int >::type print_each(print_eachSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(svilf_interal_probit(Y, H, prop, sample_adaptive, intercept, eigen_init, get_samples, tau, kappa, maxit, maxit_inner, print_each, tol));
    return rcpp_result_gen;
END_RCPP
}
// eff_lowtri
arma::vec eff_lowtri(const arma::sp_mat Y);
RcppExport SEXP _svilf_eff_lowtri(SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(eff_lowtri(Y));
    return rcpp_result_gen;
END_RCPP
}
// eff_cross_lt
arma::vec eff_cross_lt(const arma::mat Y);
RcppExport SEXP _svilf_eff_cross_lt(SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(eff_cross_lt(Y));
    return rcpp_result_gen;
END_RCPP
}
// id_lt
arma::mat id_lt(const arma::vec vv, int V);
RcppExport SEXP _svilf_id_lt(SEXP vvSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type vv(vvSEXP);
    Rcpp::traits::input_parameter< int >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(id_lt(vv, V));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_svilf_svilf_interal_logit", (DL_FUNC) &_svilf_svilf_interal_logit, 13},
    {"_svilf_svilf_interal_probit", (DL_FUNC) &_svilf_svilf_interal_probit, 13},
    {"_svilf_eff_lowtri", (DL_FUNC) &_svilf_eff_lowtri, 1},
    {"_svilf_eff_cross_lt", (DL_FUNC) &_svilf_eff_cross_lt, 1},
    {"_svilf_id_lt", (DL_FUNC) &_svilf_id_lt, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_svilf(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
