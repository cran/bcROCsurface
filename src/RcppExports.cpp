// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// asyVarVUS_C
NumericVector asyVarVUS_C(NumericVector tt, NumericMatrix D_hat, double mu_hat, NumericMatrix EstFunc, NumericMatrix Hess_inv, NumericMatrix Der_D1_hat, NumericMatrix Der_D2_hat, NumericMatrix Der_D3_hat);
RcppExport SEXP bcROCsurface_asyVarVUS_C(SEXP ttSEXP, SEXP D_hatSEXP, SEXP mu_hatSEXP, SEXP EstFuncSEXP, SEXP Hess_invSEXP, SEXP Der_D1_hatSEXP, SEXP Der_D2_hatSEXP, SEXP Der_D3_hatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tt(ttSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type D_hat(D_hatSEXP);
    Rcpp::traits::input_parameter< double >::type mu_hat(mu_hatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type EstFunc(EstFuncSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Hess_inv(Hess_invSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Der_D1_hat(Der_D1_hatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Der_D2_hat(Der_D2_hatSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Der_D3_hat(Der_D3_hatSEXP);
    rcpp_result_gen = Rcpp::wrap(asyVarVUS_C(tt, D_hat, mu_hat, EstFunc, Hess_inv, Der_D1_hat, Der_D2_hat, Der_D3_hat));
    return rcpp_result_gen;
END_RCPP
}
// vusC
double vusC(NumericVector tt, NumericMatrix dd);
RcppExport SEXP bcROCsurface_vusC(SEXP ttSEXP, SEXP ddSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tt(ttSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dd(ddSEXP);
    rcpp_result_gen = Rcpp::wrap(vusC(tt, dd));
    return rcpp_result_gen;
END_RCPP
}
