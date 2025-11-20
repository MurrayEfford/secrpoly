\name{secrpoly-package}
\alias{secrpoly-package}
\alias{secrpoly}

\docType{package}

\title{Area-search Spatially Explicit Capture--Recapture }

\description{

Functions to estimate the density and size of a spatially distributed 
animal population sampled by searching polygons or transects.
  
The package contains a subset of the functionality in \pkg{secr} that may 
in future be removed from \pkg{secr}.

}

\details{

\tabular{ll}{
Package: \tab secrpoly\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2025-11-21\cr
License: \tab GNU General Public License Version 2 or later\cr
}

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture with area searches. \emph{Ecology}
\bold{92}, 2202--2207.

Efford, M. G. (2025) The SECR book. A handbook of spatially explicit 
capture--recapture methods. Version 1.0.1. Zenodo \url{https://doi.org/10.5281/zenodo.15109938}.
Online at \url{https://murrayefford.github.io/SECRbook/}. 

}

\keyword{ package }

\seealso{

\code{\link{secrpoly.fit}}

}

\examples{

\donttest{

predict(FTHL.0)

}
}