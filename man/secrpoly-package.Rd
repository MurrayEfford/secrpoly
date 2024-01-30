\name{secrpoly-package}
\alias{secrpoly-package}
\alias{secrpoly}

\docType{package}

\title{Spatially Explicit Capture--Recapture Models }

\description{

  Functions to estimate the density and size of a spatially
  distributed animal population sampled with an array of passive
  detectors, such as traps, or by searching polygons or transects.

}

\details{

\tabular{ll}{
Package: \tab secrpoly\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2024-01-29\cr
License: \tab GNU General Public License Version 2 or later\cr
}

This contains a subset of the functionality in \pkg{secr}.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture with area searches. \emph{Ecology}
\bold{92}, 2202--2207.

}

\keyword{ package }

\seealso{

\code{\link{secrpoly.fit}}

}

\examples{

\dontrun{

FTHL.fit <- secr.fit(hornedlizardCH, buffer = 80, ncores = 7)
predict(FTHL.fit)

}
}