// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ADMM_genlasso
Rcpp::List ADMM_genlasso(Rcpp::NumericMatrix X, Rcpp::NumericVector Y, Rcpp::NumericMatrix D, double lambda, double rho, double abstol, double reltol, int maxiter);
RcppExport SEXP _augADMM_ADMM_genlasso(SEXP XSEXP, SEXP YSEXP, SEXP DSEXP, SEXP lambdaSEXP, SEXP rhoSEXP, SEXP abstolSEXP, SEXP reltolSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(ADMM_genlasso(X, Y, D, lambda, rho, abstol, reltol, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// augADMM_genlasso
List augADMM_genlasso(NumericMatrix A_mat, NumericVector b_vec, NumericMatrix D_mat, double lambda, double rho, double alpha, double abstol, double reltol, int maxiter);
RcppExport SEXP _augADMM_augADMM_genlasso(SEXP A_matSEXP, SEXP b_vecSEXP, SEXP D_matSEXP, SEXP lambdaSEXP, SEXP rhoSEXP, SEXP alphaSEXP, SEXP abstolSEXP, SEXP reltolSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A_mat(A_matSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b_vec(b_vecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type D_mat(D_matSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< double >::type reltol(reltolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(augADMM_genlasso(A_mat, b_vec, D_mat, lambda, rho, alpha, abstol, reltol, maxiter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_augADMM_ADMM_genlasso", (DL_FUNC) &_augADMM_ADMM_genlasso, 8},
    {"_augADMM_augADMM_genlasso", (DL_FUNC) &_augADMM_augADMM_genlasso, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_augADMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
