// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// gaussian_grammat_rcpp
NumericMatrix gaussian_grammat_rcpp(const NumericMatrix x, double bandwidth, int n, int d);
RcppExport SEXP _dHSIC_gaussian_grammat_rcpp(SEXP xSEXP, SEXP bandwidthSEXP, SEXP nSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussian_grammat_rcpp(x, bandwidth, n, d));
    return rcpp_result_gen;
END_RCPP
}
// discrete_grammat_rcpp
LogicalMatrix discrete_grammat_rcpp(const IntegerMatrix x, int n, int d);
RcppExport SEXP _dHSIC_discrete_grammat_rcpp(SEXP xSEXP, SEXP nSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(discrete_grammat_rcpp(x, n, d));
    return rcpp_result_gen;
END_RCPP
}
// median_bandwidth_rcpp
double median_bandwidth_rcpp(const NumericMatrix x, int n, int d);
RcppExport SEXP _dHSIC_median_bandwidth_rcpp(SEXP xSEXP, SEXP nSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(median_bandwidth_rcpp(x, n, d));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dHSIC_gaussian_grammat_rcpp", (DL_FUNC) &_dHSIC_gaussian_grammat_rcpp, 4},
    {"_dHSIC_discrete_grammat_rcpp", (DL_FUNC) &_dHSIC_discrete_grammat_rcpp, 3},
    {"_dHSIC_median_bandwidth_rcpp", (DL_FUNC) &_dHSIC_median_bandwidth_rcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_dHSIC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}