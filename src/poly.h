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

#include <R_ext/Applic.h>    // for Rdqags

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

// probability of count with distribution specified by binomN 
double countp (int count, int binomN, double lambda);
//--------------------------------------------------------------------------

double SegCircle2 (double p1x, double p1y, double p2x, double p2y, double scx, double scy, 
                   double r);

//---------------------------------------------------------------------

double expmin (double x);

double distance1 (const rpoint p1, const rpoint p2);

double rcount (const int binomN, const double lambda, const double Tsk);
//---------------------------------------------------------------------

rpoint getxy(
        const double l, 
        double cumd[], 
                   const rpoint line[], 
                                    const int kk, 
                                    const int offset);   // double before 2022-01-18

//---------------------------------------------------------------------

double randomtime (double p);
double randomtimel (double lambda);
//---------------------------------------------------------------------

void probsort (
        const int n, 
        std::vector<trap_animal> &tran);

//---------------------------------------------------------------------

double gr (
        const int fn, 
        Rcpp::NumericVector gsb, 
        const rpoint xy, 
        const rpoint animal);
//---------------------------------------------------------------------

// random point from 2-D radial distribution specified by g function 
Rcpp::NumericVector gxy (const int fn, 
                         const Rcpp::NumericVector par, 
                         const double w);

//---------------------------------------------------------------------

double hazard (double pp);

void getdetspec (
        const Rcpp::IntegerVector &detect, 
        const int fn, 
        const int nc,  
        const int nc1, 
        const int cc, 
        const int nmix, 
        const int nd, 
        const int nk, 
        const int ss, 
        const int mm, 
        const Rcpp::IntegerVector &PIA, 
        const Rcpp::NumericVector &miscparm, 
        const std::vector<int> &start, 
        std::vector<double> &detspec);

//---------------------------------------------------------------------
// 
double gpois (int count, double lambda);
double gbinom(int count, int size, double p);
double pski ( int binomN, int count, double Tski, double g, double pI);

//--------------------------------------------------------------------------

// double d2 (int k, int m, double A1[], double A2[], int A1rows, int A2rows);

double d2cpp (
        const int k, 
        const int m, 
        const Rcpp::NumericMatrix &A1, 
        const Rcpp::NumericMatrix &A2);

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

int nval(int detect0, int nc1, int cc, int ss, int nk);

void squaredistcpp (Rcpp::NumericMatrix &dist2);

bool insidecpp (
        const Rcpp::NumericVector &xy,
        const int    n1,
        const int    n2,
        const Rcpp::NumericMatrix &poly);

void fillngcpp(const int nc, 
               const int gg, 
               const Rcpp::IntegerVector &grp, 
               std::vector<int> &ng);

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

double integral1DNRcpp
    (const int fn, 
     const int m, 
     const int c, 
     const RcppParallel::RMatrix<double> &gsbval, 
     const RcppParallel::RMatrix<double> &traps,
     const RcppParallel::RMatrix<double> &mask, 
     const int n1, 
     const int n2);

//--------------------------
// original 2-D using Rdqags
//--------------------------
double integral2Dcpp  (
        const int &fn,
        const int &m,
        const int &c,
        const RcppParallel::RMatrix<double> &gsbval,
        const RcppParallel::RMatrix<double> &poly,
        const RcppParallel::RMatrix<double> &mask,
        const int &n1,
        const int &n2,
        double ex[]);

//-----------------------------------------------------
// alternative 2-D using RcppNumerical Numer::integrate
//-----------------------------------------------------
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

double hintegral1Ncpp (
        const int fn, 
        const std::vector<double> &gsb);

double hintegral2Ncpp (
        const int fn, 
        const std::vector<double> &gsb); 

#endif

