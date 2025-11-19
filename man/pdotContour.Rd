\name{pdotContour}
\alias{pdotContour}

\title{
Contour Detection Probability
}
\description{

  Display contours of the net probability of detection p.(X), or the
area within a specified distance of detectors. 

}
\usage{

pdotContour(traps, border = NULL, nx = 64, detectfn = 14,
    detectpar = list(g0 = 0.2, sigma = 25, z = 1), noccasions = NULL,
    binomN = NULL, levels = seq(0.1, 0.9, 0.1), poly =
    NULL, poly.habitat = TRUE, plt = TRUE, add = FALSE, fill = NULL, ...)

}

\arguments{
  \item{traps}{ \code{traps} object }
  \item{border}{ width of blank margin around the outermost detectors }
  \item{nx}{ dimension of interpolation grid in x-direction }
  \item{detectfn}{ integer code or character string for shape of
    detection function 14 = hazard halfnormal etc. -- see \link[secr]{detectfn}} 
  \item{detectpar}{ list of values for named parameters of detection
  function }
  \item{noccasions}{ number of sampling occasions}
  \item{binomN}{ integer code for discrete distribution (see
    \code{\link{secrpoly.fit}}) }
  \item{levels}{ vector of levels for p.(X)}
  \item{poly}{ matrix of two columns, the x and y
    coordinates of a bounding polygon (optional) }
  \item{poly.habitat}{ logical as in \code{\link[secr]{make.mask}} }
  \item{plt}{ logical to plot contours }
  \item{add}{ logical to add contour(s) to an existing plot }
  \item{fill}{ vector of colours to fill contours (optional)}
  \item{\dots}{ other arguments to pass to \code{contour} }
}
\details{

\code{pdotContour} constructs a rectangular mask and applies \code{\link[secr]{pdot}} to
compute the p.(X) at each mask point.

If \code{traps} has a \link[secr]{usage} attribute then \code{noccasions} is
  set accordingly; otherwise it must be provided.

If \code{traps} is for multiple sessions then detectpar should be a list
of the same length, one component per session, and noccasions may be a
numeric vector of the same length.

Increase \code{nx} for smoother lines, at the expense of speed.

}
\value{
  
  Coordinates of the plotted contours are returned as a list with one
  component per polygon. The list is returned invisibly if \code{plt =
    TRUE}.

  For multi-session input (\code{traps}) the value is a list of such
  lists, one per session.
  
}

\references{

  Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially
explicit capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand
\url{https://www.otago.ac.nz/density/}.

}

\note{
 
  To suppress contour labels, include the argument \code{drawlabels =
  FALSE} (this will be passed via \dots to \code{contour}). Other useful
  arguments of \code{contour} are \code{col} (colour of contour lines)
  and \code{lwd} (line width).

}

\seealso{

  \code{\link[secr]{pdot}}, \code{\link[secr]{make.mask}}
  
}
\examples{

\donttest{

FTHLtraps <- traps(hornedlizardCH)
plot(FTHLtraps, border = 100)
pdotContour(FTHLtraps, detectfn = 'HHN', nx = 128, detectpar =
    detectpar(FTHL.0), levels = c(0.1, 0.01, 0.001),
    noccasions = 14, add = TRUE)
    
}

}