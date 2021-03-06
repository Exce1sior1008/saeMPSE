// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// varner
List varner(arma::vec ni, arma::mat X, arma::vec Y, int method);
RcppExport SEXP _saeMSPE_varner(SEXP niSEXP, SEXP XSEXP, SEXP YSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ni(niSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(varner(ni, X, Y, method));
    return rcpp_result_gen;
END_RCPP
}
// mspeNERdb
arma::vec mspeNERdb(arma::vec ni, arma::mat X, arma::vec Y, arma::mat Xmean, int K, int C, int method);
RcppExport SEXP _saeMSPE_mspeNERdb(SEXP niSEXP, SEXP XSEXP, SEXP YSEXP, SEXP XmeanSEXP, SEXP KSEXP, SEXP CSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ni(niSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmean(XmeanSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type C(CSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeNERdb(ni, X, Y, Xmean, K, C, method));
    return rcpp_result_gen;
END_RCPP
}
// mspeNERjack
arma::vec mspeNERjack(arma::vec ni, arma::mat X, arma::vec Y, arma::mat Xmean, int method);
RcppExport SEXP _saeMSPE_mspeNERjack(SEXP niSEXP, SEXP XSEXP, SEXP YSEXP, SEXP XmeanSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ni(niSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmean(XmeanSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeNERjack(ni, X, Y, Xmean, method));
    return rcpp_result_gen;
END_RCPP
}
// mspeNERmcjack
arma::vec mspeNERmcjack(arma::vec ni, arma::mat X, arma::vec Y, arma::mat Xmean, int K, int method);
RcppExport SEXP _saeMSPE_mspeNERmcjack(SEXP niSEXP, SEXP XSEXP, SEXP YSEXP, SEXP XmeanSEXP, SEXP KSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ni(niSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmean(XmeanSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeNERmcjack(ni, X, Y, Xmean, K, method));
    return rcpp_result_gen;
END_RCPP
}
// mspeNERpb
arma::vec mspeNERpb(arma::vec ni, arma::mat X, arma::vec Y, arma::mat Xmean, int K);
RcppExport SEXP _saeMSPE_mspeNERpb(SEXP niSEXP, SEXP XSEXP, SEXP YSEXP, SEXP XmeanSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ni(niSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmean(XmeanSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeNERpb(ni, X, Y, Xmean, K));
    return rcpp_result_gen;
END_RCPP
}
// mspeNERsumca
arma::vec mspeNERsumca(arma::vec ni, arma::mat X, arma::vec Y, arma::mat Xmean, int K, int method);
RcppExport SEXP _saeMSPE_mspeNERsumca(SEXP niSEXP, SEXP XSEXP, SEXP YSEXP, SEXP XmeanSEXP, SEXP KSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ni(niSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmean(XmeanSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeNERsumca(ni, X, Y, Xmean, K, method));
    return rcpp_result_gen;
END_RCPP
}
// mspeFHdb
arma::vec mspeFHdb(arma::vec Y, arma::mat X, arma::vec D, int K, int C, int method);
RcppExport SEXP _saeMSPE_mspeFHdb(SEXP YSEXP, SEXP XSEXP, SEXP DSEXP, SEXP KSEXP, SEXP CSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type C(CSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeFHdb(Y, X, D, K, C, method));
    return rcpp_result_gen;
END_RCPP
}
// mspeFHjack
arma::vec mspeFHjack(arma::vec Y, arma::mat X, arma::vec D, int method);
RcppExport SEXP _saeMSPE_mspeFHjack(SEXP YSEXP, SEXP XSEXP, SEXP DSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeFHjack(Y, X, D, method));
    return rcpp_result_gen;
END_RCPP
}
// mspeFHmcjack
arma::vec mspeFHmcjack(arma::vec Y, arma::mat X, arma::vec D, int K, int method);
RcppExport SEXP _saeMSPE_mspeFHmcjack(SEXP YSEXP, SEXP XSEXP, SEXP DSEXP, SEXP KSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeFHmcjack(Y, X, D, K, method));
    return rcpp_result_gen;
END_RCPP
}
// mspeFHpb
arma::vec mspeFHpb(arma::vec Y, arma::mat X, arma::vec D, int K);
RcppExport SEXP _saeMSPE_mspeFHpb(SEXP YSEXP, SEXP XSEXP, SEXP DSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeFHpb(Y, X, D, K));
    return rcpp_result_gen;
END_RCPP
}
// mspeFHsumca
arma::vec mspeFHsumca(arma::vec Y, arma::mat X, arma::vec D, int K, int method);
RcppExport SEXP _saeMSPE_mspeFHsumca(SEXP YSEXP, SEXP XSEXP, SEXP DSEXP, SEXP KSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(mspeFHsumca(Y, X, D, K, method));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_saeMSPE_varner", (DL_FUNC) &_saeMSPE_varner, 4},
    {"_saeMSPE_mspeNERdb", (DL_FUNC) &_saeMSPE_mspeNERdb, 7},
    {"_saeMSPE_mspeNERjack", (DL_FUNC) &_saeMSPE_mspeNERjack, 5},
    {"_saeMSPE_mspeNERmcjack", (DL_FUNC) &_saeMSPE_mspeNERmcjack, 6},
    {"_saeMSPE_mspeNERpb", (DL_FUNC) &_saeMSPE_mspeNERpb, 5},
    {"_saeMSPE_mspeNERsumca", (DL_FUNC) &_saeMSPE_mspeNERsumca, 6},
    {"_saeMSPE_mspeFHdb", (DL_FUNC) &_saeMSPE_mspeFHdb, 6},
    {"_saeMSPE_mspeFHjack", (DL_FUNC) &_saeMSPE_mspeFHjack, 4},
    {"_saeMSPE_mspeFHmcjack", (DL_FUNC) &_saeMSPE_mspeFHmcjack, 5},
    {"_saeMSPE_mspeFHpb", (DL_FUNC) &_saeMSPE_mspeFHpb, 4},
    {"_saeMSPE_mspeFHsumca", (DL_FUNC) &_saeMSPE_mspeFHsumca, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_saeMSPE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
