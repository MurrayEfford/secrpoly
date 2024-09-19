###############################################################################
## package 'secrpoly'
## generalsecrloglik.R
## 2024-01-29
###############################################################################

# dettype
# 
# polygonX    = 3,
# transectX   = 4,
# polygon     = 6,
# transect    = 7,

#--------------------------------------------------------------------------------
allhistpolygon <- function (detectfn, realparval, haztemp, hk, H, pi.density, PIA, 
                           CH, xy, binomNcode, grp, usge, mask, pmixn, maskusage,
                           grain, ncores, minprob, debug=FALSE) {
  nc <- nrow(CH)
  m <- nrow(pi.density)
  s <- ncol(usge)
  nmix <- nrow(pmixn)
  ngroup <- length(levels(grp))
  sump <- numeric(nc)
  for (x in 1:nmix) {
      hx <- if (any(binomNcode==-2)) matrix(haztemp$h[x,,], nrow = m) else -1 ## lookup sum_k (hazard)
      hi <- if (any(binomNcode==-2)) haztemp$hindex else -1                   ## index to hx
      temp <- polygonhistoriescpp(
        as.integer(nc),
        as.integer(detectfn[1]),
        as.integer(grain),
        as.integer(ncores),
        as.double(minprob),          
        as.integer(binomNcode),
        as.integer(CH),   
        as.matrix(xy$xy),
        as.vector(xy$start),
        as.integer(as.numeric(grp))-1L,
        as.double(hk),
        as.double(H),
        as.matrix(realparval),
        matrix(1,nrow=s, ncol=nmix),  ## pID?
        as.matrix(mask),
        as.matrix (pi.density),
        as.integer(PIA[1,,,,x]),
        as.matrix(usge),
        as.matrix (hx),                
        as.matrix (hi),      
        as.matrix(maskusage),
        as.integer(debug)
      )
      sump <- sump + pmixn[x,] * temp
  }
  if (debug) {
    cat("SUM log(PRWI) ", sum(log(sump)), "\n")
  }
  sump
}
#--------------------------------------------------------------------------------

integralprw1poly <- function (detectfn, realparval0, haztemp, hk, H, pi.density, PIA0, 
                              CH0, binomNcode, grp, usge, mask, pmixn, maskusage,
                              grain, ncores, minprob, debug = FALSE) {
    
  nc <- dim(PIA0)[2]
  nr <- nrow(CH0)       ## unique naive animals (1 or nc)
  m <- nrow(pi.density)
  nmix <- nrow(pmixn)
  if (length(grp)<=1) grp <- rep(1,nc)
  s <- ncol(usge)
  ngroup <- length(levels(grp))
  sump <- numeric(nc)
  for (x in 1:nmix) {
      hx <- if (any(binomNcode==-2)) matrix(haztemp$h[x,,], nrow = m) else -1 ## sum_k (hazard)
      hi <- if (any(binomNcode==-2)) haztemp$hindex else -1                   ## index to hx
    for (g in 1:ngroup) {
      ok <- as.integer(grp) == g
      temp <- polygonhistoriescpp(
        as.integer(nr),
        as.integer(detectfn[1]),
        as.integer(grain),
        as.integer(ncores),
        as.double(minprob),          
        as.integer(binomNcode),
        as.integer(CH0),   
        as.matrix(0L),  # empty for null history
        as.vector(0L),  # empty for null history
        as.integer(g)-1L,
        as.double(hk),
        as.double(H),
        as.matrix(realparval0),
        matrix(1,nrow=s, ncol=nmix),  ## pID?
        as.matrix(mask),
        as.matrix (pi.density),
        as.integer(PIA0[1,1:nr,,,x]),
        as.matrix(usge),
        as.matrix (hx),                
        as.matrix (hi),      
        as.matrix(maskusage),
        as.integer(debug)
      )
      if (nr == 1) temp <- rep(temp, nc)
      sump[ok] <- sump[ok] + pmixn[x,ok] * (1-temp[ok])
    }
  }
  sump
}
#-------------------------------------------------------------------------------

################################################################################
generalsecrloglikfn <- function (
  beta, 
  parindx, 
  link, 
  fixed, 
  designD, 
  designNE, 
  design, 
  design0, 
  CL, 
  detectfn,
  learnedresponse,
  sessionlevels,
  data,
  details,
  dig = 3, betaw = 10, neglik = TRUE)
  
  # Return the negative log likelihood for spatial capture-recapture model
  
  # Transformed parameter values (density, g0, sigma, z etc.) are passed in the vector 'beta'
  # 'detectfn' is integer code for detection function
  #    0 = halfnormal, 1 = hazard, 2 = exponential etc.
  # 'CL' is logical for conditional (CL=T) vs full (CL=F) likelihood
  # details$trace=T sends a one-line report to the screen
  
{
  #--------------------------------------------------------------------------------
  sessionLL <- function (data) {
    ## log likelihood for one session
    ## in multi-session case must get session-specific data from lists
    #---------------------------------------------------
    sessnum <- data$sessnum  # changed from argument 2021-06-22
    nc1 <- max(data$nc,1)
    PIA <- design$PIA[sessnum, 1:nc1, 1:data$s, 1:data$K, ,drop = FALSE]
    PIA0 <- design0$PIA[sessnum, 1:nc1, 1:data$s, 1:data$K, ,drop = FALSE]
    ## unmodelled beta parameters, if needed
    miscparm <- getmiscparm(details$miscparm, detectfn, beta, parindx, details$cutval)
    #---------------------------------------------------
    
    density <- getmaskpar(!CL, D, data$m, sessnum, details$unmash, 
                          attr(data$capthist, 'n.mash'))
    if (CL) {
      pi.density <- matrix(1/data$m, nrow=data$m, ncol=1)  
    }
    else {
      Dsum <- apply(density,2,sum)   ## by group
      Nm <- density * getcellsize(data$mask)
      pi.density <- sweep(density, MARGIN = 2, STATS = Dsum, FUN = '/')
    }
    #---------------------------------------------------
    ## allow for scaling of detection
    Dtemp <- if (D.modelled) mean(D[,1,sessnum]) else NA
    Xrealparval <- reparameterize (realparval, detectfn, details,
                                   data$mask, data$traps, Dtemp, data$s)
    Xrealparval0 <- reparameterize (realparval0, detectfn, details,
                                    data$mask, data$traps, Dtemp, data$s)
    if (details$debug>2) browser()

    ## check valid parameter values
    if (!all(is.finite(Xrealparval))) {
      cat ('beta vector :', beta, '\n')
      warning ("extreme 'beta' in 'generalsecrloglikfn' ",
               "(try smaller stepmax in nlm Newton-Raphson?)")
      return (1e10)
    }
    ## DOES NOT ALLOW FOR GROUP VARIATION IN DENSITY
    ## more thoughts 2015-05-05
    ## could generalize by
    ## -- making Dtemp a vector of length equal rows in realparval
    ## -- matching either
    ##      first group (as before)
    ##      sum of all groups
    ##      own group [PROBLEM: locating group of each realparval row]
    ## in all cases density is the mean over mask points
    
    ## CHECK use of Dtemp in regionN.R, sim.secr.R
    ## PERHAPS for consistency make a function to construct Dtemp vector
    ## given mask, model, group matching rule (first, sum, own)
    
    #####################################################################
    pmixn <- getpmix (data$knownclass, PIA, Xrealparval)  ## membership prob by animal
    
    ## precompute gk, hk for polygon and transect detectors
    dimension <- (data$dettype[1] %in% c(3,6)) + 1   ## 1 = 1D, 2 = 2D
    # 2019-11-25 not safe to use multithreading with 2-D integration 
    # 2019-11-25 therefore using repeated 1-D integration
    convexpolygon <- is.null(details$convexpolygon) || details$convexpolygon
    gkhk <- makegkPolygoncpp (
        as.integer(detectfn), 
        as.integer(dimension), 
        as.logical(convexpolygon), 
        as.integer(details$grain), 
        as.integer(details$ncores),
        as.matrix(Xrealparval), 
        as.integer(data$cumk),
        as.matrix(data$traps), 
        as.matrix(data$mask))
    if (details$debug) {
        cat("sum(hk) ", sum(gkhk$hk), "\n")  
    }
    
    #######################################################################
    ## hazard for exclusive detectors or related
    haztemp <- gethazard (data$m, data$binomNcode, nrow(Xrealparval), gkhk$hk, PIA, data$usge)
    
    ## model detection histories (prw) conditional on detection (pdot)
    if (data$nc == 0) {
        prw <- 1  ## simple if no animals detected
    }
    else {
        prw <- allhistpolygon (detectfn, Xrealparval, haztemp, gkhk$hk, gkhk$H, pi.density, PIA, 
                               data$CH, data$xy, data$binomNcode, data$grp, data$usge, data$mask,
                               pmixn, data$maskusage, details$grain, details$ncores, details$minprob,
                               debug = details$debug>3)
    }    
    ## polygon types
    if (learnedresponse) {   ## overwrite gk,hk with model for naive animal
        gkhk <- makegkPolygoncpp (
            as.integer(detectfn), 
            as.integer(details$grain),
            as.integer(details$ncores),
            as.matrix(Xrealparval0),
            as.integer(data$cumk),
            as.matrix(data$traps),
            as.matrix(data$mask))
        if (all(data$dettype %in% c(3,4))) {
            ## hazard for exclusive detectors or related bug fix 2020-04-24
            haztemp <- gethazard (data$m, data$binomNcode, nrow(Xrealparval0), gkhk$hk, PIA0, data$usge)
        }
    }
    if (!is.null(details$externalpdot)) {
        pdot <- rep(sum(data$externalpdot * pi.density), nc1)
    }
    else {
        pdot <- integralprw1poly (detectfn, Xrealparval0, haztemp, gkhk$hk, 
                                  gkhk$H, pi.density, PIA0, data$CH0, data$binomNcode, data$grp, 
                                  data$usge, data$mask, pmixn, data$maskusage, details$grain, 
                                  details$ncores, details$minprob, debug = details$debug>3)
    }
    
    ngroup <- max(length(levels(data$grp)),1)
    comp <- matrix(0, nrow = 6, ncol = ngroup)
    for (g in 1:ngroup) {
        ok <- as.integer(data$grp) == g
        #----------------------------------------------------------------------
        
        comp[1,g] <- if (any(is.na(prw)) || any(prw<=0)) NA else sum(log(prw[ok]))
        
        #----------------------------------------------------------------------
        ## Adjust for undetected animals 
        ## or density relative.
        if (details$relativeD != 1) {
            comp[2,g] <- if (any(is.na(pdot)) || any(pdot<=0)) NA else -sum(log(pdot[ok]))
        }
        
        #----------------------------------------------------------------------
        
        if (!CL && !details$relativeD) {
            ng <- sum(ok)
            nonzero <- ng
            N <- sum(Nm[,g])
            if (ng == 0) {
                meanpdot <- pdot
            }
            else {
                meanpdot <- ng / sum(1/pdot[ok])
            }
            if (data$n.distrib == 1 && .localstuff$iter == 0 && nonzero>N) {
                warning("distribution = 'binomial' ",
                        "but number detected n (", nonzero, 
                      ") exceeds initial value of N (", round(N,1), ")")
          }
              
          comp[3,g] <- if (is.na(meanpdot) || (meanpdot <= 0)) NA 
              else switch (data$n.distrib+1,
                               dpois(nonzero, N * meanpdot, log = TRUE),
                               lnbinomial (nonzero, N, meanpdot),
                               NA)
      }
      #----------------------------------------------------------------------
      # adjustment for mixture probabilities when class known
      known <- sum(data$knownclass[ok]>1)
      if (details$nmix>1 && known>0) {
          nb <- details$nmix + 1
          nm <- tabulate(data$knownclass[ok], nbins = nb)
          pmix <- attr(pmixn, 'pmix')
          firstx <- match ((1:details$nmix)+1, data$knownclass)
          pdpmix <- pdot[firstx] * pmix
          pdpmix <- pdpmix[!is.na(pdpmix)]
          comp[4,1] <- sum(nm[-1] * log(pdpmix / sum(pdpmix)))
      }
      #----------------------------------------------------------------------
    }   ## end loop over groups
    if (details$debug>=1) {
        ## display likelihood components summed over groups, and logmultinomial constant
        comp <- apply(comp,1,sum)
        cat(comp[1], comp[2], comp[3], comp[4], comp[5], comp[6], data$logmult, '\n')
    }
    sum(comp) + data$logmult
  
  } ## end sessionLL
  
  ######################################################################################
  ## Main line of generalsecrloglikfn
  ######################################################################################
  if (details$debug>4) browser()
  nsession <- length(sessionlevels)
  #--------------------------------------------------------------------
  # Fixed beta
  beta <- fullbeta(beta, details$fixedbeta)
  #--------------------------------------------------------------------
  # Detection parameters
  detparindx <- parindx[!(names(parindx) %in% c('D', 'noneuc'))]
  detlink <- link[!(names(link) %in% c('D', 'noneuc'))]
  realparval  <- makerealparameters (design, beta, detparindx, detlink, fixed)
  realparval0 <- makerealparameters (design0, beta, detparindx, detlink, fixed)
  #--------------------------------------------------------------------
  sessmask <- lapply(data, '[[', 'mask')
  grplevels <- unique(unlist(lapply(data, function(x) levels(x$grp))))
  #---------------------------------
  # Density
  D.modelled <- !CL & is.null(fixed$D)
  if (!CL ) {
      D <- getD (designD, beta, sessmask, parindx, link, fixed,
               grplevels, sessionlevels, parameter = 'D', details$relativeD)
  }
  NE <- NULL

  loglik <- sum(sapply (data, sessionLL)) 
  .localstuff$iter <- .localstuff$iter + 1  
  if (details$trace) {
      fixedbeta <- details$fixedbeta
      if (!is.null(fixedbeta))
          beta <- beta[is.na(fixedbeta)]
      cat(format(.localstuff$iter, width=4),
          formatC(round(loglik,dig), format='f', digits=dig, width=10),
          formatC(beta, format='f', digits=dig+1, width=betaw),
          '\n')
      flush.console()
  }
  loglik <- ifelse(is.finite(loglik), loglik, -1e10)
  ifelse (neglik, -loglik, loglik)
}  ## end of generalsecrloglikfn
############################################################################################

