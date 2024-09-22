###############################################################################
## package 'secrpoly'
## pdot.R
## return net detection probability in 'traps' for home ranges centred at X
## 2024-01-29
###############################################################################

## pdot is used in --

## CVa
## CVpdot
## esa     
## fxTotal (fxTotal.secrpoly uses secrpoly pdot)
## make.mask (pdotmin option restricted to point detectors)
## reparameterize.esa (restricted to point detectors)
## bias.D  (restricted to point detectors)
## pdot.contour


pdot <- function (X, traps, detectfn = 14, detectpar = list(g0 = 0.2, sigma = 25, z = 1),
                  noccasions = NULL, binomN = NULL, userdist = NULL, ncores = NULL) {

    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coord in col 2

    ncores <- setNumThreads(ncores)
    grain <- if (ncores==1) 0 else 1

    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if ((detectfn > 9) && (detectfn<14) && is.null(detectpar$cutval))
        stop ("requires 'cutval' for detectfn 10:13")
    if (ms(traps))
        stop ("requires single-session traps")

    truncate <- ifelse(is.null(detectpar$truncate), 1e+10, detectpar$truncate)

    detectpars <- unlist(detectpar[parnames(detectfn)])
    if ((detectfn>9) && (detectfn<14))  detectpars <- c(detectpars, detectpar$cutval)
    if (length(detectpars)<3) detectpars <- c(detectpars,0)
    miscparm <- numeric(4);   ## dummy

    if (!is.null(usage(traps))) {
        usge <- usage(traps)
        if (is.null(noccasions)) {
            noccasions <- ncol(usage(traps))
        }
        else {
            if (noccasions < ncol(usage(traps))) {
                warning ("specified noccasions less than ncol of usage matrix")
            }
            if (noccasions > ncol(usage(traps)))
                stop ("specified noccasions exceeds ncol of usage matrix")
        }
    }
    else {
        if (is.null(noccasions))
            stop("must specify noccasions when traps does not have usage attribute")
        usge <- matrix(1, ndetector(traps), noccasions)
    }
    dettype <- detectorcode(traps, noccasions = noccasions)
    binomN <- getbinomN (binomN, detector(traps))
    markocc <- markocc(traps)
    
    if (is.null(markocc)) markocc <- rep(1,noccasions)
    if (!inherits(X, 'mask')) {
        X <- matrix(unlist(X), ncol = 2)
    }
    if (all(detector(traps) %in% .localstuff$polydetectors)) {
        if (!is.null(userdist))
            stop("userdist incompatible with polygon-like detectors")
        if (!(detectfn %in% 14:20))
            stop("pdot requires hazard detectfn for polygon-type detectors")
        k <- table(polyID(traps))   ## also serves transectID
        K <- length(k)              ## number of polygons/transects
        k <-  c(k,0)                ## zero terminate for 
        cumk <- cumsum(c(0,table(polyID(traps))))
        convexpolygon <- TRUE
        dim <- if (any(detector(traps) %in% c('transect', 'transectX'))) 1 else 2
            
        warning("assuming convex polygons in pdot()")
        temp <- hdotpolycpp (
          as.matrix(X),
          as.matrix(traps),
          as.matrix(usge),
          as.integer(markocc),
          as.integer(cumk),
          as.integer(detectfn),
          as.double(detectpars),
          as.logical(convexpolygon),
          as.integer(dim),
          as.integer(grain),
          as.integer(ncores))
        1 - exp(-temp)   ## probability detected at least once, given total hazard
    }
    else stop ("expecting polygon detectors")
    
}
############################################################################################

pdot.contour <- function (traps, border = NULL, nx = 64, detectfn = 14,
                          detectpar = list(g0 = 0.2, sigma = 25, z = 1),
                          noccasions = NULL, binomN = NULL,
                          levels = seq(0.1, 0.9, 0.1),
                          poly = NULL, poly.habitat = TRUE, plt = TRUE, add = FALSE, fill = NULL, ...) {

    if (ms(traps)) {
        if (length(noccasions) == 1)
            noccasions <- rep(noccasions,length(traps))
        output <- mapply(pdot.contour, traps, detectpar=detectpar, noccasions=noccasions,
                         MoreArgs = list(border = border, nx = nx,
                         detectfn = detectfn, binomN = binomN,
                         levels = levels, poly = poly, poly.habitat = poly.habitat, plt = plt, add = add, ...))
        if (plt)
            invisible(output)
        else
            output
    }
    else {
        if (is.null(border))
            border <- 5 * spatialscale(detectpar, detectfn)
        tempmask <- make.mask (traps, border, nx = nx, type = 'traprect')
        xlevels <- unique(tempmask$x)
        ylevels <- unique(tempmask$y)
        binomN <- getbinomN (binomN, detector(traps))
        z <- pdot(tempmask, traps, detectfn, detectpar, noccasions, binomN)
        if (!is.null(poly)) {
            OK <- pointsInPolygon(tempmask, poly)
            if (poly.habitat)
                z[!OK] <- 0
            else
                z[OK] <- 0
        }
        if (plt) {
            contour (xlevels, ylevels, matrix(z, nrow = nx), add = add, levels = levels, ...)


            ## optional fillin 2015-05-15
            if (!is.null(fill)) {
                z[z < (0.999 * min(levels))] <- NA
                levels <- c(0,levels,1)
                .filled.contour(xlevels, ylevels,  matrix(z, nrow = nx), levels= levels,
                                col = fill)
            }


            invisible(contourLines(xlevels, ylevels, matrix(z, nrow = nx), levels = levels))
        }
        else
            contourLines(xlevels, ylevels, matrix(z, nrow = nx), levels = levels)
    }
}
############################################################################################

