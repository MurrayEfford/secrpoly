############################################################################################
## package 'secrpoly'
## esa.R
## 2024-09-20
############################################################################################

esa.secrpoly <- function (object, sessnum = 1, beta = NULL, real = NULL, 
                          noccasions = NULL, ncores = NULL, Dweight = FALSE, ...)

# Return vector of 'a' for given g0, sigma, [z (if hazard fn) ] and session
# detectfn is integer code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential
# 'sessnum' is integer index of session (factor level of the 'session' attribute in capthist)
# object must have at least capthist, mask, detectfn

## strictly doesn't need data, so better to avoid when object not available...
{
    if (ms(object))
        capthists <- object$capthist[[sessnum]]
    else
        capthists <- object$capthist

    if (nrow(capthists) == 0) {
        return (numeric(0))    
    }

    # Remainder applies only to polydetectors
    if (ms(object$mask)) {
        mask <- object$mask[[sessnum]]
    }
    else {
        mask <- object$mask
    }

    if (is.null(beta) & is.null(real)) {
        beta <- object$fit$par
    }
    
    fullbeta <- fullbeta(beta, object$details$fixedbeta)
    trps     <- traps(capthists)  ## need session-specific traps
    n        <- max(nrow(capthists), 1)
    s        <- ncol(capthists)
    dettype  <- secr:::secr_detectorcode(trps, noccasions = s)
    constant <- !is.null(noccasions)    ## fix 2011-04-07
    if (is.null(noccasions)) {
        noccasions <- s
    }
    markocc <- markocc(traps(capthists))
    if (is.null(markocc))
        markocc <- rep(1,s)  ## simple marking
    allsighting <- !any(markocc>0)
    if (allsighting) {
        ## drop all zero histories, consider sighting as if marking
        # capthists <- subset(capthists, apply(capthists!=0,1,sum)>0)
        warning("use of allsighting here untested")
        ## markocc[] <- 1
    }
    else {
        #not right, but getting there
        #object$capthist <- subset(object$capthist, occasions = (markocc>0))
    }
    #----------------------------------------------------------------------
    nmix    <- getnmix (object$details)
    knownclass <- getknownclass(capthists, nmix, object$hcov)
    k <- getk(trps)  # number of vertices per detector, zero terminated
    K <- if (length(k)>1) length(k)-1 else k  # number of detectors
    m      <- length(mask$x)            ## need session-specific mask...
    cellsize <- getcellsize(mask)       ## length or area
    
    if (is.null(object$model$D) || !Dweight) {
        D <- rep(1,m)
    }
    else {
        if (!is.null(beta)) object$fit$par <- beta
        D <- secr:::secr_predictD(object, mask, group = NULL, session = sessnum, 
                           parameter = 'D')
    }
    pi.density <- matrix(D/sum(D), ncol = 1)       # dim m x 1
    
    # Non-Euclidean distance parameters
    NElist <- secr:::secr_makeNElist(object, object$mask, group = NULL, sessnum)
  
    #----------------------------------------------------------------------
    if (constant) {
        ## assume constant
        if (is.null(beta))
            realparval0 <- detectpar(object)
        else {
            realparval0 <- secr:::secr_makerealparameters (object$design0, fullbeta,
                object$parindx, object$link, object$fixed)  # naive
            realparval0 <- as.list(realparval0[1,])
            realparval0$cutval <- attr(object$capthist,'cutval')  ## 2016-05-22 may be NULL
        }
        a <- cellsize * sum(pdot(X = mask, traps = trps, detectfn = object$detectfn,
                             detectpar = realparval0, noccasions = noccasions, 
                             ncores = ncores))
        out <- rep(a,n)
    }
    else {
        if (is.null(beta)) {
            if (is.null(real))
                stop ("requires real parameter values")
            PIA0 <- array(1, dim=c(1,n,s,K,nmix))
            realparval0 <- matrix(rep(real, rep(n,length(real))), nrow = n)   ## UNTRANSFORMED
        }
        else {
            if (object$CL) {
                PIA0 <- object$design0$PIA[sessnum,1:n,1:s,1:K,,drop = FALSE]
            }
            else {
                ## use whatever columns are present as the number of groups
                ## should be constant over sessions
                PIA0 <- object$design0$PIA[sessnum,,1:s,1:K,,drop = FALSE]
                ## fill array with PI appropriate to grouping of i-th animal
                PIA0 <- PIA0[1,group.factor(capthists, object$groups),,,,drop = FALSE]
            }

            realparval0 <- makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive

        }
        ## not compatible with sigmak parameterizations
        Dtemp <- NA
        Xrealparval0 <- secr:::secr_reparameterize (realparval0, object$detectfn, object$details,
                                        mask, trps, Dtemp, s)

        usge <- usage(trps)
        if (is.null(usge) || object$details$ignoreusage) {
                usge <- matrix(1, nrow = K, ncol = s)
            used <- 1
        }
        else {
            ## force to binary 2012-12-17
            used <- (usge > 1e-10) * 1
        }
        if (any(used == 0))
            PIA0 <- PIA0 * rep(rep(t(used),rep(n,s*K)),nmix)

        miscparm <- numeric(4)

        ## not if first detector type is telemetry
        if (dettype[1] %in% c(13)) {
            out <- NA
        }
        else {
          ncores <- setNumThreads(ncores)  
          grain <- if (ncores==1) 0 else 1
          gkhk <- makegk (dettype, object$detectfn, trps, mask, object$details, sessnum, NElist, 
                          pi.density, miscparm, Xrealparval0, grain, ncores)
          if (any(dettype==0)) {
              CH0 <- secr:::secr_nullCH (c(n,s), object$design0$individual)
          }
          else {
              CH0 <- secr:::secr_nullCH (c(n,s,K), object$design0$individual)
          }
          binomNcode <- recodebinomN(dettype, object$details$binomN, telemcode(trps))
          pmixn <- getpmix (knownclass, PIA0, Xrealparval0)
          pdot <- integralprw1poly (
              detectfn    = object$detectfn,
              realparval0 = Xrealparval0, 
              haztemp     = secr:::secr_gethazard(m, binomNcode, nrow(Xrealparval0), gkhk$hk, PIA0, usge), 
              hk          = gkhk$hk, 
              H           = gkhk$H, 
              pi.density  = pi.density, 
              PIA0        = PIA0, 
              CH0         = CH0, 
              binomNcode  = binomNcode, 
              grp         = rep(1,n),
              usge        = usge, 
              mask        = mask, 
              pmixn       = pmixn, 
              grain       = grain, 
              ncores      = ncores,
              minprob     = object$details$minprob)
          
          out <- pdot * cellsize * sum(D)
        }
    }
    out
}
############################################################################################
