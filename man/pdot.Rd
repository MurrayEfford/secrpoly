\name{pdotPoly}
\alias{pdotPoly}

\title{ Net Detection Probability }

\description{

Compute spatially explicit net probability of detection for individual(s) 
at given coordinates.

}

\usage{
pdotPoly(X, traps, detectfn = 14, detectpar = list(g0 = 0.2,
    sigma = 25, z = 1), noccasions = NULL, binomN = NULL,
    userdist = NULL, ncores = NULL) 

}

\arguments{
  \item{X}{ vector or 2-column matrix of coordinates }
  \item{traps}{ \code{traps} object }
  \item{detectfn}{ integer code for detection function q.v. }
  \item{detectpar}{ a named list giving a value for each parameter of detection function }
  \item{noccasions}{ number of sampling intervals (occasions) }
  \item{binomN}{ integer code for discrete distribution (see
    \code{\link{secrpoly.fit}}) }
  \item{userdist}{user-defined distance function or matrix (see \link[secr]{userdist})}
  \item{ncores}{integer number of threads}
}
\details{

  If \code{traps} has a \link[secr]{usage} attribute then \code{noccasions} is
  set accordingly; otherwise it must be provided.
  
The probability computed is \eqn{p.(\mathbf{X}) = 1 - \prod\limits _{k}
\{1 - p_s(\mathbf{X},k)\}^{S}}{p.(X) = 1 - (1 - prod(p_s(X,k))^S} where
the product is over the detectors in \code{traps}, excluding any not
used on a particular occasion. The per-occasion detection function
\eqn{p_s}{p_s} is halfnormal (0) by default, and is assumed not to vary
over the \eqn{S} occasions.

From 4.6.11, the detection parameters g0, lambda0 and sigma for point detectors 
may be detector- and occasion-specific. This is achieved by providing a vector 
of values that is replicated internally to fill a matrix with dimensions 
ntraps x noccasions (i.e. in trap order for occasion 1, then occasion 2 etc.)

For detection functions (10) and (11) the signal threshold `cutval' should be
included in \code{detectpar}, e.g., \code{detectpar = list(beta0 = 103, beta1
= -0.11, sdS = 2, cutval = 52.5)}.

The calculation is not valid for single-catch traps because
\eqn{p.(\mathbf{X})}{p.(X)} is reduced by competition between animals.

\code{userdist} cannot be set if `traps' is any of polygon, polygonX,
transect or transectX. if \code{userdist} is a function requiring
covariates or values of parameters `D' or `noneuc' then \code{X} must
have a covariates attribute with the required columns.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link[secr]{setNumThreads}}).

}

\value{

A vector of probabilities, one for each row in X.

}

\seealso{ 
  \code{\link[secr]{make.mask}},
  \code{\link[secr]{Detection functions}}, 
  \code{\link{pdotPolyContour}}
}

\examples{

\donttest{
  temppoly <- make.poly()
  ## per-session detection probability for a centred individual.
  ## Using default noccasions = 5.
  pdotPoly (c(0,0), temppoly, detectpar = list(g0 = 0.2, sigma = 25),
    noccasions = 5)
}

}
