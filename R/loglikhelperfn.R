###############################################################################
## package 'secrpoly'
## loglikhelperfn.R
## 2024-01-29
###############################################################################

getmiscparm <- function(miscparm, detectfn, beta, parindx, cutval) {
    ## miscparm is used to package beta parameters that are not modelled
    ## and hence do not have a beta index specified by parindx.
    ## This includes the signal threshold and the mean and sd of noise.
    ## miscparm is passed to the C++ likelihood code, and 
    ## also as mask attribute to userdistfn
    nmiscparm <- length(miscparm)
    miscparm <- numeric(max(4, nmiscparm)) 
    if (detectfn %in% c(12,13))            ## experimental signal-noise
        miscparm[1:3] <- c(cutval, beta[max(unlist(parindx))+(1:2)])   ## fudge: last 2
    else if (detectfn %in% c(10,11))        ## Dawson & Efford 2009 models
        miscparm[1] <- cutval
    else if (nmiscparm > 0)
        miscparm[1:nmiscparm] <- beta[max(unlist(parindx)) + (1:nmiscparm)]
    miscparm
}
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
getD <- function (designD, beta, mask, parindx, link, fixed,
                  grouplevels, sessionlevels, parameter = 'D',
                  relativeD = FALSE) {
    ## apply to either 'D' or 'noneuc'
    if (!is.function(designD)) {
        if ((is.null(designD) || nrow(designD)==0) && (is.null(fixed[[parameter]]))) return(NULL)
    }
    if (ms(mask))
        nmask <- max(sapply(mask, nrow))
    else
        nmask <- nrow(mask)
    ngroup <- max(1, length(grouplevels))
    nsession <- length(sessionlevels)
    D <- array(dim = c(nmask, ngroup, nsession))
    dimnames(D) <- list(1:nrow(D), grouplevels, sessionlevels)
    if (!is.null(fixed[[parameter]])) {
        D[,,] <- fixed[[parameter]]
    }
    else {
        if (is.function(designD)) {
            ## fixed 2021-12-10, 2022-05-24
            if (ms(mask)) {
                for (session in 1:nsession) {
                    m <- nrow(mask[[session]])
                    D[1:m,,session] <- designD(beta[parindx[[parameter]]], mask[[session]], ngroup, 1)
                }
            } 
            else {
                m <- nrow(mask)
                D[1:m,,1] <- designD(beta[parindx[[parameter]]], mask, ngroup, 1)
            }
        }
        else {
            beta <- beta[parindx[[parameter]]]
            Dfn <- attr(designD, 'Dfn')
            if (is.function(Dfn)) {
                    D[,,] <- Dfn(designD, beta)
            }
            else {
                D[,,] <- designD %*% beta  # linear predictor
            }
            D[,,] <- untransform (D, link[[parameter]])
        }
        # silently truncate D at zero
        # allow non-positive noneuc
        if (parameter == 'D')
            D[D<0] <- 0
    }
    D
}
#--------------------------------------------------------------------------------
getmaskpar <- function(OK, D, m, sessnum, unmash, nmash) {
    if (!OK) {
        NULL   ## not in model
    }
    else {
        sessnum <- min(dim(D)[3],sessnum)
        density <- matrix(D[1:m,,sessnum], nrow = m)
        if (!all(is.finite(density))) {
            cat ('densities :', head(density), '\n')
            ## 2022-03-20
            warning ("bad densities in 'secrloglikfn' ",
                "(try different optimisation method, link, or model?)")
            # stop ("bad densities in 'secrloglikfn' ",
            #     "(try different optimisation method, link, or model?)")
        }
        ## optional scaling by session-specific number of clusters
        if (unmash) {
            if (!is.null(nmash))
                density <- density * length(nmash)
        }
        density
    }
}
#--------------------------------------------------------------------------------
makegk <- function(dettype, detectfn, trps, mask, details, sessnum, noneuc, 
                   D, miscparm, realparval, grain, ncores) {
    ## precompute gk, hk for polygon and transect detectors
    ## k-1 because we have zero-terminated these vectors
    ## noneuc not used 2025-11-20
    k <- getk(trps)
    K <- if (length(k)>1) length(k)-1 else k
    cumk <- cumsum(c(0,k))[1:length(k)]
    dimension <- (dettype[1] %in% c(3,6)) + 1   ## 1 = 1D, 2 = 2D
    convexpolygon <- is.null(details$convexpolygon) || details$convexpolygon
    gkhk <- makegkPolygoncpp (
        as.integer(detectfn), 
        as.integer(dimension), 
        as.logical(convexpolygon), 
        as.integer(grain), 
        as.integer(ncores), 
        as.matrix(realparval), 
        as.integer(cumk),
        as.matrix(trps), 
        as.matrix(mask))
    
    gkhk
}
