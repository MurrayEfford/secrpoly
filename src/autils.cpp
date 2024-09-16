#include "poly.h"

using namespace std;
using namespace Rcpp;

double minimumexp = -100;

//--------------------------------------------------------------------------

double expmin (double x)
{
    if (x < minimumexp)
        return(0);
    else
        return(exp(x));
}
//--------------------------------------------------------------------------

// index to vector element corresponding to cell i,j,k in 3D array
// stored in column-major order 

int i3 (int i, int j, int k, int ii, int jj) {
    return(ii * (jj * k + j) + i);
}
//--------------------------------------------------------------------------

// index to vector element corresponding to cell i,j,k,l in 4D array
// stored in column-major order 

int i4 (int i, int j, int k, int l, int ii, int jj, int kk) {
    return (ii *(jj*(kk*l + k) + j) + i);
}
//--------------------------------------------------------------------------

// customised dbinom 
double gbinom(int count, int size, double p)
{
    double x, q;
    int i;
    if ((count < 0) || (count > 0 && p <= 0)) {
        x = 0;
    }
    else if (count == 0) {
        q = 1 - p;
        x = q;
        for (i=1; i< size; i++) x *= q;
    }
    else {
        boost::math::binomial_distribution<> bin(size, p);
        x = boost::math::pdf(bin, count);
    }
    return (x);   
}
//--------------------------------------------------------------------------

// use distance1 rather than distance because of clash with std::distance
double distance1 (const rpoint p1, const rpoint p2) {
    return(std::sqrt ((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y)));
}
//--------------------------------------------------------------------------

double zrcpp (double r, int detectfn, NumericVector par)
{
    if (detectfn == 14) {  // hazard halfnormal
        return (exp(-r*r / 2 / par(1) / par(1)));    
    }
    else {
        if (detectfn == 15) {  // hazard hazard rate
            return (1 - exp(- pow(r /par(1), - par(2))));
        }
        else if (detectfn == 16) {  // hazard exponential
            return (exp(-r / par(1)));
        }
        else if (detectfn == 17) {  // hazard annular normal
            return (exp(-(r-par(2))*(r-par(2)) / 
                    2 / par(1)/ par(1)));
        }
        else if (detectfn == 18) {  // hazard cumulative gamma
            // return (R::pgamma(r,par(2),par(1)/par(2),0,0)); 
            boost::math::gamma_distribution<> gam(par(2),par(1)/par(2));
            return (boost::math::cdf(complement(gam,r))); 
        }
        else if (detectfn == 19) {  // hazard variable power
            return (exp(- pow(r /par(1), par(2))));
        }
        else 
            return (R_NaN);  //Rcpp::stop("unknown or invalid detection function in zrcpp"));
    }
}

// [[Rcpp::export]]
bool insidecpp (
        const NumericVector &xy,
        const int    n1,
        const int    n2,
        const NumericMatrix &poly)
{
    // Is point xy inside poly?
    // Based on contribution on s-news list by Peter Perkins 23/7/96
    // We assume poly is closed, and in col-major order (x's then y's)
    
    double theta = 0;
    double cutoff = 1e-6;
    int k;
    int ns;
    double N;
    double d;
    ns = n2 - n1 + 1;   // number of selected points 
    std::vector<double> temp((ns+1) * 2);
    
    // get & translate to coords centered at each test point 
    for (k=0; k < ns; k++)
    {
        temp[k]      = poly(k + n1,0) - xy[0];    // x 
        temp[k + ns] = poly(k + n1,1) - xy[1];    // y 
    }
    
    for (k=0; k < (ns-1); k++)
    {
        N = temp[k] * temp[k+1 + ns] - temp[k + ns] * temp[k+1];
        d = temp[k] * temp[k+1]      + temp[k + ns] * temp[k+1 + ns];
        if (fabs(d)>0) { N = N/fabs(d);  d = d/fabs(d); }
        theta += atan2(N, d);
    }
    theta = fabs(theta);
    return (fabs(theta - 2* M_PI) < cutoff);    // M_PI is cmath.h constant 
}
//--------------------------------------------------------------------------


// Calculate the length of intersection of a line segment and a circle
// Based on C code of Paul Bourke November 1992
// Line segment is defined from p1 to p2
// Circle is of radius r and centred at sc
// Two potential points of intersection given by
// p = p1 + mu1 (p2-p1)
// p = p1 + mu2 (p2-p1)
// Return 0 if line segment does not intersect circle

double SegCircle2 (
        double p1x, double p1y, 
        double p2x, double p2y, 
        double scx, double scy, 
        double r
) 
{
    double a,b,c;
    double bb4ac;
    double dpx;
    double dpy;
    
    double mu1;
    double mu2;
    int p1in;
    int p2in;
    
    double i1x;
    double i1y;
    double i2x;
    double i2y;
    
    int i1between;
    double d1,d2;
    double seg = 0;
    
    // case where both p1 and p2 inside circle 
    
    // Rprintf ("p1 %6.3f %6.3f\n", p1x, p1y);
    // Rprintf ("p2 %6.3f %6.3f\n", p2x, p2y);
    // Rprintf ("sc %6.3f %6.3f\n", scx, scy);
    // Rprintf ("r %6.3f \n", r);
    
    p1in = ((scx - p1x) * (scx - p1x) + 
        (scy - p1y) * (scy - p1y)) < (r * r);
    p2in = ((scx - p2x) * (scx - p2x) + 
        (scy - p2y) * (scy - p2y)) < (r * r);
    if (p1in && p2in) {        
        seg = std::sqrt ((p1x - p2x) * (p1x - p2x) + 
            (p1y - p2y) * (p1y - p2y));
        return (seg);
    }
    
    dpx = p2x - p1x;
    dpy = p2y - p1y;
    
    a = dpx * dpx + dpy * dpy;
    b = 2 * (dpx * (p1x - scx) + dpy * (p1y - scy));
    c = scx * scx + scy * scy;
    c += p1x * p1x + p1y * p1y;
    c -= 2 * (scx * p1x + scy * p1y);
    c -= r * r;
    bb4ac = b * b - 4 * a * c;
    
    // case of no intersection 
    if ((fabs(a) < 1e-10) || (bb4ac < 0)) {
        return (0);   
    }
    
    mu1 = (-b + std::sqrt(bb4ac)) / (2 * a);
    mu2 = (-b - std::sqrt(bb4ac)) / (2 * a);
    
    i1x = p1x + mu1 * (p2x - p1x);
    i1y = p1y + mu1 * (p2y - p1y);
    i2x = p1x + mu2 * (p2x - p1x);
    i2y = p1y + mu2 * (p2y - p1y);
    
    if (((mu1<0) && (mu2<0)) || ((mu1>1) && (mu2>1))) {
        // no intersection 
        seg = 0;
    }
    else {
        if (((mu1<0) && (mu2>1)) || ((mu1>1) && (mu2<0))) {
            // both inside 
            seg = std::sqrt ((p1x - p2x) * (p1x - p2x) + 
                (p1y - p2y) * (p1y - p2y));
        }
        else {
            if ((mu1>0) && (mu1<1) && (mu2>0) && (mu2<1)) {
                // two intersections 
                seg = std::sqrt ((i1x - i2x) * (i1x - i2x) + 
                    (i1y - i2y) * (i1y - i2y));
            }
            else {
                // one intersection 
                d1 = std::sqrt((i1x - p1x) * (i1x * p1x) + 
                    (i1y - p1y) * (i1y - p1y));
                d2 = std::sqrt((i1x - p2x) * (i1x * p2x) + 
                    (i1y - p2y) * (i1y - p2y));
                i1between = std::sqrt(a) < (d1 + d2 + 1e-10);
                if (p1in) {
                    if (i1between) {
                        i2x = p1x;
                        i2y = p1y;
                    }
                    else {
                        i1x = p1x;
                        i1y = p1y;
                    }
                }
                if (p2in) {
                    if (i1between) {
                        i2x = p2x;
                        i2y = p2y;
                    }
                    else {
                        i1x = p2x;
                        i1y = p2y;
                    }
                }
                seg = std::sqrt ((i1x - i2x) * (i1x - i2x) + 
                    (i1y - i2y) * (i1y - i2y));
            }
        }
    }
    return(seg);    
}

//----------------------------------------------------------------

// return probability g(r) for given detection function fn 
// used in simsecr.cpp and trapping.cpp 
// double gr (
//         const int fn,
//         const Rcpp::NumericVector gsb,
//         const rpoint xy,
//         const rpoint animal) {
//     double r;
//     fnptrC fnp;
//     fnp = getgfns(fn);
//     r = distance1 (xy, animal);
//     return (fnp(as<std::vector<double>>(gsb),r));
// }
//----------------------------------------------------------------


// double hazard (double pp) {
//     if (pp > (1-fuzz))  // pp close to 1.0 - approx limit 
//         pp = huge;      // g0 very large (effecti inf hazard) 
//     else {
//         if (pp <= 0) 
//             pp = 0;
//         else 
//             pp = -log(1-pp);
//     }
//     return(pp);
// }
//=============================================================

// detect may take values -
// 0  multi-catch traps
// 1  binary proximity detectors
// 2  count  proximity detectors
// 3  exclusive polygon detector
// 4  exclusive transect detector
// 5  signal detector
// 6  polygon detector
// 7  transect detector
// 8  times  (undocumented)
// 9  cue    (undocumented) -- removed in secr 2.10.0
// 12 signalnoise


//--------------------------------------------------------------------------

// Functions to characterize detector type 
// polygon, transect and signal detector types must be constant 
// across occasions 

bool anyexclusive (const IntegerVector detect) {
    bool exclusive = false;
    for (int s=0; s< detect.size(); s++) {
        if ((detect[s]==0) || (detect[s]==3) || (detect[s]==4))
            exclusive = true;
    }
    return exclusive;
}

bool anycapped  (const IntegerVector detect) {
    bool capped = false;
    for (int s=0; s<detect.size(); s++) {
        if (detect[s]==8)
            capped = true;
    }
    return capped;
}

bool anypolygon  (const IntegerVector detect) {
    bool polygon = false;
    for (int s=0; s<detect.size(); s++) {
        if ((detect[s]==3) || (detect[s]==6) )
            polygon = true;
    }
    return polygon;
}

bool anytransect (const IntegerVector detect) {
    bool transect = false;
    for (int s=0; s<detect.size(); s++) {
        if ((detect[s]==4) || (detect[s]==7))
            transect = true;
    }
    return transect;
}

bool anysignal (const IntegerVector detect) {
    bool signal = false;
    for (int s=0; s<detect.size(); s++) {
        if ((detect[s]==5) || (detect[s]==12))
            signal = true;
    }
    return signal;
}

bool anytelemetry (const IntegerVector detect) {
    bool telemetry = false;
    for (int s=0; s<detect.size(); s++) {
        if (detect[s]==13)
            telemetry = true;
    }
    return telemetry;
}

//  check if we need to consider variation among individuals 
// i.e. check if detection parameters constant for given s,k 
bool anyvarying (
        const int    nc,     // number of capture histories (or groups if PIA0 has that dim) 
        const int    ss,     // number of occasions 
        const int    nk,     // number of traps 
        const int    nmix,   // number of mixture classes 
        const IntegerVector &PIA0  // lookup which g0/sigma/b combination to use for given n, S, K [naive] 
) {
    int i,n,s,k,x;
    int wxi;
    bool indiv = false;
    for (s=0; s<ss; s++) {
        for (k=0; k<nk; k++) {
            for (x=0; x<nmix; x++) {
                wxi = i4(0,s,k,x,nc,ss,nk);       
                i = PIA0[wxi];
                for (n=1; n<nc; n++) {
                    wxi = i4(n,s,k,x,nc,ss,nk);    
                    if (i != PIA0[wxi]) {
                        indiv = true; break;
                    }
                }
            }
        }
    }
    return(indiv);
}
//--------------------------------------------------------------------

bool alltelemetry (const IntegerVector detect) {
    bool telemetry = true;
    for (int s=0; s<detect.size(); s++) {
        if ((detect[s]!=13))
            telemetry = false;
    }
    return telemetry;
}

bool allpoint (const IntegerVector detect, bool allowsignal, bool allowtelem) {
    bool point;
    bool OK = true;
    for (int s=0; s<detect.size(); s++) {
        point = (detect[s]==0) || (detect[s]==1) || (detect[s]==2) || detect[s] == 8
        || (detect[s]==10) || (detect[s]==11)
        || (allowsignal && ((detect[s]==5) || (detect[s]==12)))
        || (allowtelem && ((detect[s]==13)));
        OK = OK && point;
    }
    return OK;
}

bool allcapped  (const IntegerVector detect) {
    bool OK = true;
    for (int s=0; s<detect.size(); s++) {
        OK = OK && (detect[s] == 8);
    }
    return OK;
}

bool allmulti (const IntegerVector detect) {
    bool notmulti = false;
    for (int s=0; s<detect.size(); s++) {
        if (detect[s]!=0)
            notmulti = true;
    }
    return (!notmulti);
}
//--------------------------------------------------------------------

// Do parameter values for naive animals differ at all from those for other animals ?
bool anyb (const NumericMatrix &gsbval, const NumericMatrix &gsb0val) {
    bool identical = true;
    for (int i=0; i<gsbval.size(); i++) {
        if (gsbval[i] != gsb0val[i]) identical = false;
    }
    return (!identical);
}
//==============================================================================


// probability of count for session s, detector k, animal i
// The argument 'g' is understood to be a cumulative hazard if binomN=0,
// a probability otherwise

double pski ( int binomN,
              int count,
              double Tski,
              double g,
              double pI) {
    
    double lambda;
    double result = 1.0;
    
    if (binomN == -1) {                              // binary proximity detectors : Bernoulli
        if (abs(Tski-1) > 1e-10) {                   // effort not unity; adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        if (count>0)                                 
            result = g*pI;  
        else 
            result = 1 - g*pI;
    }
    else if (binomN == 0) {                          // count detectors : Poisson 
        lambda = Tski * g * pI;
        if ((count < 0) || (count>0 && lambda<=0)) {         
            result = 0;
        }
        else if (count == 0) {
            result = exp(-lambda);            // routinely apply Tsk adjustment to cum. hazard 
        }
        else {
            boost::math::poisson_distribution<> pois(lambda);
            result = boost::math::pdf(pois,count);
        }
    }
    else if (binomN == 1) {                          // count detectors : Binomial, size from Tsk
        result = gbinom (count, round(Tski), g*pI); 
    }
    else if (binomN > 1) {                           // count detectors : Binomial, specified size 
        if (abs(Tski-1) > 1e-10) {                   // effort not unity, adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        result = gbinom (count, binomN, g*pI);
    }
    else result = NAN; // Rcpp::stop("binomN < -1 not allowed");  // code multi -2 separately
    
    return (result);
}
//--------------------------------------------------------------------------

double classmembership (
        const int n, 
        const int x, 
        const IntegerVector &knownclass, 
        const std::vector<double> &pmixn, 
        const int nmix) {
    
    // Return probability individual n belongs to class x. This may be binary 
    //   (0/1) in the case of known class, or continuous if class is unknown 
    
    double pmixnx = 0;
    int knownx = -1;
    
    if (knownclass[n] == 1) 
        knownx = -1;                         // unknown class 
    else
        knownx = R::imax2(0, knownclass[n]-2);  // known class 
    
    // unknown class : weighted by probability of membership  
    // known class does not require probability of membership 
    if (knownx < 0)
        pmixnx = pmixn[nmix * n + x];
    else if (knownx == x)
        pmixnx = 1.0;
    else 
        pmixnx = 0.0;
    return (pmixnx);
    
}
//--------------------------------------------------------------------------

