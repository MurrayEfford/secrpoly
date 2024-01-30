// #include <Rcpp.h>
#include "secr.h"
using namespace Rcpp;

//     Detection functions
//     2019-07-29 C++ and consolidate
//     2019-10-13 eliminate unused functions
//     2024-01-29 hazard only for polygons
// 
// The function pfnS() selects gfnr on the fly.
// pfnS() is used in naiveRPSV, simdetect, 
// 
// All g--- functions return the detection probability g(x). For count detectors this
// must be transformed to the cumulative hazard -log(1-g(x)) (see function 'hazard' 
// in utils.c).
// 
// The gfns is used 
// 
// (i) for polygon and transect detectors, 
// (ii) for trappingXXX simulations  (via pfnS)
// 
// for which the function must be calculated on the fly within the integration algorithm 
// rather than from a lookup precomputed for a fixed set of trap-mask distances.
// 
// For each form there is a corresponding selection function:
//  getzfnr   polygon.cpp
//  getgfns   pdot.cpp, trapping.cpp (pfnS)
//
//                                    getzfnr    getgfns
//
//                                    fnptr      fnptrC  
//                                              
// fn 14 hazard halfnormal            zhhnr     ghhns      
// fn 15 hazard hazard rate           zhhrr     ghhrs      
// fn 16 hazard exponential           zhexr     ghexs      
// fn 17 hazard annular normal        zhanr     ghans      
// fn 18 hazard cumulative gamma      zhcgr     ghcgs      
// fn 19 hazard variable power        etc.

//--------------------------------------------------------------------
// define functions with third parameter z

int par3 (int fn) {
    if ((fn == 15) || (fn==17) || (fn == 18))
	return(1);
    else
	return(0);
}

// hazard halfnormal 
double zhhnr (const NumericVector& param, const double r) {
    return(param[0] * exp(- r * r / 2 / param[1] / param[1]));
}
//--------------------------------------------------------------------

// hazard hazard rate 
double zhhrr (const NumericVector& param, const double r) {
    return(param[0] * ( 1 - exp(- pow(r / param[1], -param[2]))));
}
//--------------------------------------------------------------------

// hazard exponential 
double zhexr (const NumericVector& param, const double r) {
    return (param[0] * exp(-r / param[1]));
}
//--------------------------------------------------------------------

// hazard annular normal 
double zhanr (const NumericVector& param, const double r) {
    return (param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
				   param[1] / param[1]));
}
//--------------------------------------------------------------------

// hazard cumulative gamma 
double zhcgr (const NumericVector& param, const double r) {
    return ((1 - exp( - param[0] * exp(-r / param[1]))));
}

// hazard variable power 
double zhvpr (const NumericVector& param, const double r) {
    return (param[0] * exp( - pow(r / param[1], param[2]))) ;
}
//--------------------------------------------------------------------

// hazard halfnormal 
double zhhnrC (const std::vector<double>& param, const double r) {
    return(param[0] * exp(- r * r / 2 / param[1] / param[1]));
}
//--------------------------------------------------------------------

// hazard hazard rate 
double zhhrrC (const std::vector<double>& param, const double r) {
    return(param[0] * ( 1 - exp(- pow(r / param[1], -param[2]))));
}
//--------------------------------------------------------------------

// hazard exponential 
double zhexrC (const std::vector<double>& param, const double r) {
    return (param[0] * exp(-r / param[1]));
}
//--------------------------------------------------------------------

// hazard annular normal 
double zhanrC (const std::vector<double>& param, const double r) {
    return (param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
            param[1] / param[1]));
}
//--------------------------------------------------------------------

// hazard cumulative gamma 
double zhcgrC (const std::vector<double>& param, const double r) {
    return ((1 - exp( - param[0] * exp(-r / param[1]))));
}

// hazard variable power 
double zhvprC (const std::vector<double>& param, const double r) {
    return (param[0] * exp( - pow(r / param[1], param[2]))) ;
}


//====================================================================

// hazard halfnormal 
double ghhns (const std::vector<double>& param, const double r) {
    // Rprintf("%6.3f %6.3f %6.3f \n", r, param[0], param[1]); 
    return(1 - exp( - param[0] * exp(- r * r / 2 / param[1] / param[1])));
}
//--------------------------------------------------------------------

// hazard hazard rate 
double ghhrs (const std::vector<double>& param, const double r) {
    return(1 - exp( - param[0] * ( 1 - exp(- pow(r / param[1], -param[2])))));
}
//--------------------------------------------------------------------

// hazard exponential 
double ghexs (const std::vector<double>& param, const double r) {
    return (1 - exp( - param[0] * exp(-r / param[1])));
}
//--------------------------------------------------------------------

// hazard annular normal 
double ghans (const std::vector<double>& param, const double r) {
    return (1 - exp( - param[0] * exp(-(r-param[2])*(r-param[2]) / 2 /
            param[1] / param[1])));
}
//--------------------------------------------------------------------

// hazard cumulative gamma 
double ghcgs (const std::vector<double>& param, const double r) {
    return (1 - exp( - (1 - exp( - param[0] * exp(-r / param[1])))));
}
//--------------------------------------------------------------------

// hazard variable power 
double ghvps (const std::vector<double>& param, const double r) {
    return(1 - exp( - param[0] * exp(- pow(r / param[1], param[2]))));
}

//====================================================================

double pfnS (
        const int fn,
        const double d2val,
        const std::vector<double> &gsb,
        const std::vector<double> &miscparm,
        const double w2)
{
    double p = -1;
    fnptrC gfns;
    std::vector<double> tmp(4);
    
    if (d2val > w2) 
        p = 0;
    else {
        
        gfns = getgfns (fn);
        tmp[0] = gsb[0];
        tmp[1] = gsb[1];
        tmp[2] = gsb[2];
        tmp[3] = miscparm[0];
        p = gfns (tmp, std::sqrt(d2val));
    }
    return (p);
}
//--------------------------------------------------------------------

fnptrC getgfns (const int fn) 
{
    if (fn == 14)
        return(ghhns);
    else if (fn == 15)
        return(ghhrs);
    else if (fn == 16)
        return(ghexs);
    else if (fn == 17)
        return(ghans);
    else if (fn == 18)
        return(ghcgs);
    else if (fn == 19)
        return(ghvps);
    else // Rcpp::stop("unknown or invalid detection function");
        return(ghhns);
}
//--------------------------------------------------------------------

// For polygon integral 
fnptr getzfnr (const int fn) 
{
  if (fn == 14)
    return(zhhnr);
  else if (fn == 15)
    return(zhhrr);
  else if (fn == 16)
    return(zhexr);
  else if (fn == 17)
    return(zhanr);
  else if (fn == 18)
    return(zhcgr);
  else if (fn == 19)
      return(zhvpr);
  else // Rcpp::stop("unknown or invalid detection function");
  return(zhhnr);
}
//--------------------------------------------------------------------

// For polygon integral 
fnptrC getzfnrC (const int fn) 
{
    if (fn == 14)
        return(zhhnrC);
    else if (fn == 15)
        return(zhhrrC);
    else if (fn == 16)
        return(zhexrC);
    else if (fn == 17)
        return(zhanrC);
    else if (fn == 18)
        return(zhcgrC);
    else if (fn == 19)
        return(zhvprC);
    else // Rcpp::stop("unknown or invalid detection function");
    return(zhhnrC);
}
//--------------------------------------------------------------------
