// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// edist2cpp
NumericMatrix edist2cpp(const NumericMatrix& A1, const NumericMatrix& A2);
RcppExport SEXP _secrpoly_edist2cpp(SEXP A1SEXP, SEXP A2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A2(A2SEXP);
    rcpp_result_gen = Rcpp::wrap(edist2cpp(A1, A2));
    return rcpp_result_gen;
END_RCPP
}
// xydist2cpp
NumericMatrix xydist2cpp(const NumericMatrix& A1, const NumericMatrix& A2);
RcppExport SEXP _secrpoly_xydist2cpp(SEXP A1SEXP, SEXP A2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type A2(A2SEXP);
    rcpp_result_gen = Rcpp::wrap(xydist2cpp(A1, A2));
    return rcpp_result_gen;
END_RCPP
}
// nearestcpp
List nearestcpp(const NumericMatrix& xy, const NumericMatrix& traps, bool non_zero);
RcppExport SEXP _secrpoly_nearestcpp(SEXP xySEXP, SEXP trapsSEXP, SEXP non_zeroSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type xy(xySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type traps(trapsSEXP);
    Rcpp::traits::input_parameter< bool >::type non_zero(non_zeroSEXP);
    rcpp_result_gen = Rcpp::wrap(nearestcpp(xy, traps, non_zero));
    return rcpp_result_gen;
END_RCPP
}
// insidecpp
bool insidecpp(const NumericVector& xy, const int n1, const int n2, const NumericMatrix& poly);
RcppExport SEXP _secrpoly_insidecpp(SEXP xySEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP polySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type xy(xySEXP);
    Rcpp::traits::input_parameter< const int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type poly(polySEXP);
    rcpp_result_gen = Rcpp::wrap(insidecpp(xy, n1, n2, poly));
    return rcpp_result_gen;
END_RCPP
}
// getdenomcpp
Rcpp::List getdenomcpp(int fn, Rcpp::NumericVector miscparm, Rcpp::NumericMatrix mask, int mm, double sigma, double z);
RcppExport SEXP _secrpoly_getdenomcpp(SEXP fnSEXP, SEXP miscparmSEXP, SEXP maskSEXP, SEXP mmSEXP, SEXP sigmaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type fn(fnSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type miscparm(miscparmSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(getdenomcpp(fn, miscparm, mask, mm, sigma, z));
    return rcpp_result_gen;
END_RCPP
}
// gethcpp
List gethcpp(int nc1, int cc, int nmix, int nk, int ss, int mm, const IntegerVector PIA, const NumericMatrix Tsk, const NumericVector hk);
RcppExport SEXP _secrpoly_gethcpp(SEXP nc1SEXP, SEXP ccSEXP, SEXP nmixSEXP, SEXP nkSEXP, SEXP ssSEXP, SEXP mmSEXP, SEXP PIASEXP, SEXP TskSEXP, SEXP hkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nc1(nc1SEXP);
    Rcpp::traits::input_parameter< int >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type nmix(nmixSEXP);
    Rcpp::traits::input_parameter< int >::type nk(nkSEXP);
    Rcpp::traits::input_parameter< int >::type ss(ssSEXP);
    Rcpp::traits::input_parameter< int >::type mm(mmSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA(PIASEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type hk(hkSEXP);
    rcpp_result_gen = Rcpp::wrap(gethcpp(nc1, cc, nmix, nk, ss, mm, PIA, Tsk, hk));
    return rcpp_result_gen;
END_RCPP
}
// makegkPolygoncpp
List makegkPolygoncpp(const int detectfn, const int dim, const bool convex, const int grain, const int ncores, const NumericMatrix& gsbval, const IntegerVector& cumk, const NumericMatrix& traps, const NumericMatrix& mask);
RcppExport SEXP _secrpoly_makegkPolygoncpp(SEXP detectfnSEXP, SEXP dimSEXP, SEXP convexSEXP, SEXP grainSEXP, SEXP ncoresSEXP, SEXP gsbvalSEXP, SEXP cumkSEXP, SEXP trapsSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type detectfn(detectfnSEXP);
    Rcpp::traits::input_parameter< const int >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< const bool >::type convex(convexSEXP);
    Rcpp::traits::input_parameter< const int >::type grain(grainSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type gsbval(gsbvalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cumk(cumkSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type traps(trapsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(makegkPolygoncpp(detectfn, dim, convex, grain, ncores, gsbval, cumk, traps, mask));
    return rcpp_result_gen;
END_RCPP
}
// makelookupcpp
List makelookupcpp(const NumericMatrix& x);
RcppExport SEXP _secrpoly_makelookupcpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(makelookupcpp(x));
    return rcpp_result_gen;
END_RCPP
}
// hdotpolycpp
NumericVector hdotpolycpp(const NumericMatrix& xy, const NumericMatrix& traps, const NumericMatrix& Tsk, const IntegerVector& markocc, const IntegerVector& cumk, const int& detectfn, const NumericVector& gsb, const bool& convex, const int& dim, const int& grain, const int& ncores);
RcppExport SEXP _secrpoly_hdotpolycpp(SEXP xySEXP, SEXP trapsSEXP, SEXP TskSEXP, SEXP markoccSEXP, SEXP cumkSEXP, SEXP detectfnSEXP, SEXP gsbSEXP, SEXP convexSEXP, SEXP dimSEXP, SEXP grainSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type xy(xySEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type traps(trapsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type markocc(markoccSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type cumk(cumkSEXP);
    Rcpp::traits::input_parameter< const int& >::type detectfn(detectfnSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type gsb(gsbSEXP);
    Rcpp::traits::input_parameter< const bool& >::type convex(convexSEXP);
    Rcpp::traits::input_parameter< const int& >::type dim(dimSEXP);
    Rcpp::traits::input_parameter< const int& >::type grain(grainSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(hdotpolycpp(xy, traps, Tsk, markocc, cumk, detectfn, gsb, convex, dim, grain, ncores));
    return rcpp_result_gen;
END_RCPP
}
// ontransectcpp
bool ontransectcpp(NumericVector xy, NumericMatrix transect, int n1, int n2, double tol);
RcppExport SEXP _secrpoly_ontransectcpp(SEXP xySEXP, SEXP transectSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xy(xySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type transect(transectSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(ontransectcpp(xy, transect, n1, n2, tol));
    return rcpp_result_gen;
END_RCPP
}
// alongtransectcpp
double alongtransectcpp(NumericVector xy, NumericMatrix transect, int n1, int n2, double tol);
RcppExport SEXP _secrpoly_alongtransectcpp(SEXP xySEXP, SEXP transectSEXP, SEXP n1SEXP, SEXP n2SEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xy(xySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type transect(transectSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(alongtransectcpp(xy, transect, n1, n2, tol));
    return rcpp_result_gen;
END_RCPP
}
// polygonhistoriescpp
NumericVector polygonhistoriescpp(const int nc, const int detectfn, const int grain, const int ncores, const double minp, const IntegerVector binomN, const IntegerVector w, const NumericMatrix xy, const IntegerVector start, const IntegerVector group, const NumericVector hk, const NumericVector H, const NumericMatrix gsbval, const NumericMatrix pID, const NumericMatrix mask, const NumericMatrix density, const IntegerVector PIA, const NumericMatrix Tsk, const NumericMatrix h, const IntegerMatrix hindex, const LogicalMatrix mbool, const int debug);
RcppExport SEXP _secrpoly_polygonhistoriescpp(SEXP ncSEXP, SEXP detectfnSEXP, SEXP grainSEXP, SEXP ncoresSEXP, SEXP minpSEXP, SEXP binomNSEXP, SEXP wSEXP, SEXP xySEXP, SEXP startSEXP, SEXP groupSEXP, SEXP hkSEXP, SEXP HSEXP, SEXP gsbvalSEXP, SEXP pIDSEXP, SEXP maskSEXP, SEXP densitySEXP, SEXP PIASEXP, SEXP TskSEXP, SEXP hSEXP, SEXP hindexSEXP, SEXP mboolSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const int >::type detectfn(detectfnSEXP);
    Rcpp::traits::input_parameter< const int >::type grain(grainSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const double >::type minp(minpSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type binomN(binomNSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type hk(hkSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type H(HSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type gsbval(gsbvalSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type pID(pIDSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type density(densitySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA(PIASEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type h(hSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type hindex(hindexSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type mbool(mboolSEXP);
    Rcpp::traits::input_parameter< const int >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(polygonhistoriescpp(nc, detectfn, grain, ncores, minp, binomN, w, xy, start, group, hk, H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, mbool, debug));
    return rcpp_result_gen;
END_RCPP
}
// polygonfxicpp
NumericVector polygonfxicpp(const int nc, const int detectfn, const int grain, const int ncores, const double minp, const IntegerVector binomN, const IntegerVector w, const NumericMatrix xy, const IntegerVector start, const IntegerVector group, const NumericVector hk, const NumericVector H, const NumericMatrix gsbval, const NumericMatrix pID, const NumericMatrix mask, const NumericMatrix density, const IntegerVector PIA, const NumericMatrix Tsk, const NumericMatrix h, const IntegerMatrix hindex, const LogicalMatrix mbool);
RcppExport SEXP _secrpoly_polygonfxicpp(SEXP ncSEXP, SEXP detectfnSEXP, SEXP grainSEXP, SEXP ncoresSEXP, SEXP minpSEXP, SEXP binomNSEXP, SEXP wSEXP, SEXP xySEXP, SEXP startSEXP, SEXP groupSEXP, SEXP hkSEXP, SEXP HSEXP, SEXP gsbvalSEXP, SEXP pIDSEXP, SEXP maskSEXP, SEXP densitySEXP, SEXP PIASEXP, SEXP TskSEXP, SEXP hSEXP, SEXP hindexSEXP, SEXP mboolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const int >::type detectfn(detectfnSEXP);
    Rcpp::traits::input_parameter< const int >::type grain(grainSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    Rcpp::traits::input_parameter< const double >::type minp(minpSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type binomN(binomNSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type hk(hkSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type H(HSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type gsbval(gsbvalSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type pID(pIDSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type density(densitySEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type PIA(PIASEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type Tsk(TskSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type h(hSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type hindex(hindexSEXP);
    Rcpp::traits::input_parameter< const LogicalMatrix >::type mbool(mboolSEXP);
    rcpp_result_gen = Rcpp::wrap(polygonfxicpp(nc, detectfn, grain, ncores, minp, binomN, w, xy, start, group, hk, H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, mbool));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_secrpoly_edist2cpp", (DL_FUNC) &_secrpoly_edist2cpp, 2},
    {"_secrpoly_xydist2cpp", (DL_FUNC) &_secrpoly_xydist2cpp, 2},
    {"_secrpoly_nearestcpp", (DL_FUNC) &_secrpoly_nearestcpp, 3},
    {"_secrpoly_insidecpp", (DL_FUNC) &_secrpoly_insidecpp, 4},
    {"_secrpoly_getdenomcpp", (DL_FUNC) &_secrpoly_getdenomcpp, 6},
    {"_secrpoly_gethcpp", (DL_FUNC) &_secrpoly_gethcpp, 9},
    {"_secrpoly_makegkPolygoncpp", (DL_FUNC) &_secrpoly_makegkPolygoncpp, 9},
    {"_secrpoly_makelookupcpp", (DL_FUNC) &_secrpoly_makelookupcpp, 1},
    {"_secrpoly_hdotpolycpp", (DL_FUNC) &_secrpoly_hdotpolycpp, 11},
    {"_secrpoly_ontransectcpp", (DL_FUNC) &_secrpoly_ontransectcpp, 5},
    {"_secrpoly_alongtransectcpp", (DL_FUNC) &_secrpoly_alongtransectcpp, 5},
    {"_secrpoly_polygonhistoriescpp", (DL_FUNC) &_secrpoly_polygonhistoriescpp, 22},
    {"_secrpoly_polygonfxicpp", (DL_FUNC) &_secrpoly_polygonfxicpp, 21},
    {NULL, NULL, 0}
};

RcppExport void R_init_secrpoly(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
