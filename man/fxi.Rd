\name{fxi}
\alias{fxi.secrpoly}
\alias{fxi}

\title{ Probability Density of Home Range Centre }

\description{
  Compute the probability density function for the estimated
  location of one or more range centres at particular points X. The pdf is given by
  \eqn{f(X_j|\omega_i) = \mbox{Pr}(\omega_i|X_j)\pi(X_j)}{f(X_j|\omega_i) =
  Pr(\omega_i|X)\pi(X_j)}, where \eqn{\pi(X)} is the probability density
  of range centres across the mask (Borchers and Efford 2008).
  
  This function is a method for objects of class `secrpoly'.
  `fxi' is a generic exported by \pkg{secr}. See \code{\link{fxiContour}}, 
  \code{\link{fxiMode}} and \code{\link{fxTotal}} for applications.
  
}

\usage{

\method{fxi}{secrpoly}(object, i = NULL, sessnum = 1, X = NULL, ncores = NULL, ...)

}

\arguments{
  \item{object}{ a fitted secr model }
  \item{i}{ integer or character vector of individuals (defaults to all in fxi.secr),
  or a single individual as input to fxi.mode }
  \item{sessnum}{ session number if \code{object$capthist} spans
  multiple sessions}
  \item{X}{2-column matrix of x- and y- coordinates (defaults to mask)}
  \item{ncores}{ integer number of threadss to be used for parallel processing}
  \item{\dots}{additional arguments (not used)}
}

\details{

  \code{fxi} computes the probability density for one or more
  detection histories; \code{X} may contain coordinates for one or
  several points; a dataframe or vector (x then y) will be coerced to a
  matrix.

  See \code{\link{fxTotal}} for a surface summed across individuals.

  If \code{i} is character it will be matched to row names of
  object$capthist (restricted to the relevant session in the case of a
  multi-session fit); otherwise it will be interpreted as a row number.
  
  Values of the pdf are normalised by dividing by the
  integral of \eqn{\mbox{Pr}(\omega_i|X)\pi(X)}{Pr(\omega_i|X)\pi(X)}
  over the habitat mask in \code{object}. (May differ in secr 4.0).
  
  Setting \code{ncores = NULL} uses the existing value from the environment variable 
  RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

}

\value{
  
Vector of probability densities

}
\references{

  Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

}

\seealso{

\code{\link{pdot.contour}}, \code{\link[graphics]{contour}}, \code{\link{fxiContour}}, 
  \code{\link{fxiMode}}, \code{\link{fxTotal}}

}

\examples{

\donttest{

plot(hornedlizardCH, border = 20)
points(550,520, pch = 16) # pick a point
fxi(FTHL.fit, i = 1, X = c(550,520))

## contour first 5 detection histories
fxiContour (FTHL.fit, i = 1:5, add = TRUE,
    plotmode = TRUE, drawlabels = FALSE)

## extract modes only
## these are slower but more reliable than those from fit.mode called directly 
## as they use a contour-based approximation for the starting point

fxiout <- fxiContour (FTHL.fit, i = 1:5, plt = FALSE, fitmode = TRUE)
t(sapply(fxiout, "[[", "mode"))

## using fill colours
## lty = 0 suppresses contour lines
## nx = 256 ensures smooth outline
plot(traps(hornedlizardCH), border = 20)
fxiContour(FTHL.fit, i = 1:5, add = TRUE, p = c(0.5,0.95), drawlabels
    = FALSE, nx = 256, fill = topo.colors(4), lty = 0)

## output as simple features
sf <- fxiContour(FTHL.fit, i = 1:3, plt = FALSE, p = c(0.5,0.95),
    nx = 256, output = 'sf', fitmode = TRUE)

## save as ESRI shapefile testsf.shp etc.
library(sf)
st_write(sf, 'testsf.shp')
## plot contours and modes
plot(st_as_sfc(sf))    # outline only
points(sf$modex, sf$modey)

## output as SpatialPolygonsDataFrame
spdf <- fxiContour(FTHL.fit, i = 1:3, plt = FALSE, p = c(0.5,0.95),
    nx = 256, output = 'SPDF', fitmode = TRUE)
sp::plot(spdf)
points(data.frame(spdf))

}

}
