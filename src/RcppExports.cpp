// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rrefBal
void rrefBal(NumericMatrix& M);
RcppExport SEXP _IneqCube_rrefBal(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type M(MSEXP);
    rrefBal(M);
    return R_NilValue;
END_RCPP
}
// onestepfastflightcube
NumericVector onestepfastflightcube(NumericVector prob, NumericMatrix Bm);
RcppExport SEXP _IneqCube_onestepfastflightcube(SEXP probSEXP, SEXP BmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Bm(BmSEXP);
    rcpp_result_gen = Rcpp::wrap(onestepfastflightcube(prob, Bm));
    return rcpp_result_gen;
END_RCPP
}
// flightphase
NumericVector flightphase(NumericVector prob, NumericMatrix Xbal);
RcppExport SEXP _IneqCube_flightphase(SEXP probSEXP, SEXP XbalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xbal(XbalSEXP);
    rcpp_result_gen = Rcpp::wrap(flightphase(prob, Xbal));
    return rcpp_result_gen;
END_RCPP
}
// flightphase_arma
arma::vec flightphase_arma(arma::mat X, arma::vec pik, double EPS);
RcppExport SEXP _IneqCube_flightphase_arma(SEXP XSEXP, SEXP pikSEXP, SEXP EPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< double >::type EPS(EPSSEXP);
    rcpp_result_gen = Rcpp::wrap(flightphase_arma(X, pik, EPS));
    return rcpp_result_gen;
END_RCPP
}
// onestepflightphase_arma
arma::vec onestepflightphase_arma(arma::mat B, arma::vec pik, double EPS);
RcppExport SEXP _IneqCube_onestepflightphase_arma(SEXP BSEXP, SEXP pikSEXP, SEXP EPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< double >::type EPS(EPSSEXP);
    rcpp_result_gen = Rcpp::wrap(onestepflightphase_arma(B, pik, EPS));
    return rcpp_result_gen;
END_RCPP
}
// rref
void rref(arma::mat& M);
RcppExport SEXP _IneqCube_rref(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type M(MSEXP);
    rref(M);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_IneqCube_rrefBal", (DL_FUNC) &_IneqCube_rrefBal, 1},
    {"_IneqCube_onestepfastflightcube", (DL_FUNC) &_IneqCube_onestepfastflightcube, 2},
    {"_IneqCube_flightphase", (DL_FUNC) &_IneqCube_flightphase, 2},
    {"_IneqCube_flightphase_arma", (DL_FUNC) &_IneqCube_flightphase_arma, 3},
    {"_IneqCube_onestepflightphase_arma", (DL_FUNC) &_IneqCube_onestepflightphase_arma, 3},
    {"_IneqCube_rref", (DL_FUNC) &_IneqCube_rref, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_IneqCube(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
