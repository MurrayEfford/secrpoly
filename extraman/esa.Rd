\name{esa}
\alias{esa}

\title{ Effective Sampling Area }

\description{
Compute derived `effective sampling area' of spatially explicit capture-recapture model (in the sense of Borchers and Efford 2008).
}

\usage{

esa(object, sessnum = 1, beta = NULL, real = NULL, noccasions = NULL, 
    ncores = NULL)

}

\arguments{
  \item{object}{ \code{secrpoly} object output from \code{secrpoly.fit}, or an
  object of class \code{c("secrlist", "list")} }
  \item{sessnum}{ index of session in object$capthist for which output required }
  \item{beta}{ vector of fitted parameters on transformed (link) scale }
  \item{real}{ vector of `real' parameters } 
  \item{noccasions}{ integer number of sampling occasions (see Details) } 
  \item{ncores}{integer }
}

\details{

\code{esa} is used by \code{derived} to compute individual-specific
effective sampling areas: \deqn{a_i (\hat{\theta}) = \int _A \:
p.(\mathbf{X};\mathbf{z}_i, \mathbf{\hat{\theta}}) \; \mathrm{d}
\mathbf{X}}{a_i = integral p.(X; z_i, theta_i) dX} where
\eqn{p.(\mathbf{X})}{p.(X)} is the probability an individual at X is
detected at least once and the \eqn{\mathbf{z}_i}{z_i} are optional
individual covariates. Integration is over the area \eqn{A} of the
habitat mask.

The argument \code{noccasions} may be used to vary the number of
sampling occasions; it works only when detection parameters are constant
across individuals and across time.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

The vector of detection parameters for \code{esa} may be specified via
\code{beta} or \code{real}, with the former taking precedence. If
neither is provided then the fitted values in \code{object$fit$par} are
used. Specifying \code{real} parameter values bypasses the various
linear predictors. Strictly, the `real' parameters are for a naive
capture (animal not detected previously).

}

\value{

Vector of animal-specific esa.

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum likelihood methods for capture--recapture studies. \emph{Biometrics}, \bold{64}, 377--385.

}

\seealso{

  \code{\link{predict.secr}},
  \code{\link{print.secr}},
  \code{\link{secrpoly.fit}},

}

\examples{

\donttest{

## fitted lambda0, sigma
esa(FTHL.0)

## force different lambda0, sigma
esa(FTHL.0, real = c(0.2, 25))
#Error in haztemp$h[x, , ] : subscript out of bounds
}

}
