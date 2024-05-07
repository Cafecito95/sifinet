// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// EstNull
Rcpp::List EstNull(arma::vec x, double gamma);
RcppExport SEXP _SiFINeT_EstNull(SEXP xSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(EstNull(x, gamma));
    return rcpp_result_gen;
END_RCPP
}
// cal_coexp
arma::mat cal_coexp(arma::mat X, arma::mat X_subcohort);
RcppExport SEXP _SiFINeT_cal_coexp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_coexp(X, X_subcohort));
    return rcpp_result_gen;
END_RCPP
}
// cal_coexp_sp
arma::mat cal_coexp_sp(arma::sp_mat X, arma::sp_mat X_subcohort);
RcppExport SEXP _SiFINeT_cal_coexp_sp(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_coexp_sp(X, X_subcohort));
    return rcpp_result_gen;
END_RCPP
}
// cal_conn
List cal_conn(arma::mat data, double thres, int m, bool abso, int niter);
RcppExport SEXP _SiFINeT_cal_conn(SEXP dataSEXP, SEXP thresSEXP, SEXP mSEXP, SEXP absoSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< bool >::type abso(absoSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_conn(data, thres, m, abso, niter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SiFINeT_EstNull", (DL_FUNC) &_SiFINeT_EstNull, 2},
    {"_SiFINeT_cal_coexp", (DL_FUNC) &_SiFINeT_cal_coexp, 1},
    {"_SiFINeT_cal_coexp_sp", (DL_FUNC) &_SiFINeT_cal_coexp_sp, 1},
    {"_SiFINeT_cal_conn", (DL_FUNC) &_SiFINeT_cal_conn, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_SiFINeT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
