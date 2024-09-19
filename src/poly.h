// include guard
#ifndef __poly_h_INCLUDED__   // if poly.h hasn't been included yet...
#define __poly_h_INCLUDED__   // #define this so the compiler knows it has been included

// see https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/stat_tut/weg/error_eg.html
// and https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/pol_tutorial/changing_policy_defaults.html
// return NAN for invalid inputs
#define BOOST_MATH_DOMAIN_ERROR_POLICY ignore_error
#include <boost/math/distributions.hpp>       // gamma, normal, lognormal distributions

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// next two lines must be in order (RcppNumerical precedes other)
#include <RcppNumerical.h>

#include <Rcpp.h>
#include <RcppParallel.h>

//==============================================================================
    // from secr

using namespace Rcpp;
using namespace RcppParallel;

// constants
#define fuzz 1e-200
#define huge 1e10
#define maxnpoly 1000   
#define maxnmix 2    
#define maxvertices 400

//-------------------
    // data structures   
//-------------------
    struct trap_animal {
        int     trap;
        int     animal;
        double  time;
    };
struct rpoint {
    double x;
    double y;
};

//--------------------------------------------------------------------------
    
    int i3 (int i, int j, int k, int ii, int jj);
int i4 (int i, int j, int k, int l, int ii, int jj, int kk);

//------------------------------------------------------
    // detectfn.cpp 
//------------------------------------------------------
    typedef double (*fnptr)(const Rcpp::NumericVector&, const double);
typedef double (*fnptrC)(const std::vector<double>&, const double);
fnptr getzfnr (int fn);
fnptrC getgfns (int fn);
fnptrC getzfnrC (int fn);

double pfnS (
    const int fn,
    const double d2val,
    const std::vector<double> &gsb,
    const std::vector<double> &miscparm,
    const double w2);

//--------------------------------------------------------------------------
    
    double SegCircle2 (double p1x, double p1y, double p2x, double p2y, double scx, double scy, 
                       double r);

//---------------------------------------------------------------------
    
    double expmin (double x);

double distance1 (const rpoint p1, const rpoint p2);

// double gr (
    //         const int fn, 
    //         Rcpp::NumericVector gsb, 
    //         const rpoint xy, 
    //         const rpoint animal);
//---------------------------------------------------------------------
    
    // double hazard (double pp);

double gbinom(int count, int size, double p);
double pski ( int binomN, int count, double Tski, double g, double pI);

//--------------------------------------------------------------------------
    
    // Functions to characterize detector type 
// polygon, transect and signal detector types must be constant across occasions

bool anyexclusive (const Rcpp::IntegerVector detect);
bool anycapped (const Rcpp::IntegerVector detect);
bool anypolygon (const Rcpp::IntegerVector detect);
bool anytransect (const Rcpp::IntegerVector detect);
bool anysignal (const Rcpp::IntegerVector detect);
bool anytelemetry (const Rcpp::IntegerVector detect);
bool alltelemetry (const Rcpp::IntegerVector detect);
bool allpobool (const Rcpp::IntegerVector detect, bool allowsignal, bool allowtelem);
bool allcapped  (const Rcpp::IntegerVector detect);
bool allmulti (const Rcpp::IntegerVector detect);
bool allpoint (const Rcpp::IntegerVector detect, bool allowsignal, bool allowtelem);

bool anyvarying (const int nc, const int ss, const int nk, const int nmix,
                 const Rcpp::IntegerVector &PIA0);
bool anyb (
    const Rcpp::NumericMatrix &gsbval, 
    const Rcpp::NumericMatrix &gsb0val);

// miscellaneous functions

bool insidecpp (
    const Rcpp::NumericVector &xy,
    const int    n1,
    const int    n2,
    const Rcpp::NumericMatrix &poly);

//---------------------------------------------------------------
    // Return probability individual n belongs to class x. This may be binary 
//   (0/1) in the case of known class, or continuous if class is unknown 
double classmembership (
    const int n, 
    const int x, 
    const Rcpp::IntegerVector &knownclass, 
    const std::vector<double> &pmixn, 
    const int nmix);

void yab(double x[], int *i, int *np, double poly[], double *a, double *b);
void fy(double *x, int n, void *ex);
void fx(double *x, int n, void *ex);
void fx1 (double *x, int n, void *ex);

//==============================================================================
    
    
    rpoint getxycpp(
        const double l, 
        const std::vector<double> &cumd, 
        const RcppParallel::RMatrix<double> &line, 
        const int n1, 
        const int n2);

//-------------------------------------
// see 'secr' for original using Rdqags 
//-------------------------------------
    
//-----------------------------------------------------
// alternative 2-D using RcppNumerical Numer::integrate
//-----------------------------------------------------
    
    double hintegral1DNRcpp (
        const int fn, 
        const std::vector<double> &gsb);

double hintegral2DNRcpp (
    const int fn, 
    const std::vector<double> &gsb); 

double integral1DNRcpp (
    const int fn, 
    const int m, 
    const int c, 
    const RcppParallel::RMatrix<double> &gsbval, 
    const RcppParallel::RMatrix<double> &traps,
    const RcppParallel::RMatrix<double> &mask, 
    const int n1, 
    const int n2);

double integral2DNRcpp  (
    const int &fn,
    const int &m,
    const int &c,
    const RcppParallel::RMatrix<double> &gsbval,
    const RcppParallel::RMatrix<double> &poly,
    const RcppParallel::RMatrix<double> &mask,
    const int &n1,
    const int &n2,
    const bool &convex);

#endif