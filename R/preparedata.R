###############################################################################
## package 'secrpoly'
## preparedata.R
## 2024-01-29
###############################################################################

#--------------------------------------------------------------------------------
getk <- function(traps) {
  if (!all(detector(traps) %in% .localstuff$polydetectors)) {
    nrow(traps)
  }
  else {
    if (all(detector(traps) %in% c('polygon','polygonX'))) {        
      k <- table(polyID(traps))       
    }
    else if (all(detector(traps) %in% c('transect','transectX'))) {        
      k <- table(transectID(traps))   # transectX, transect
    }
    else stop ("unrecognised poly detector type")
    c(k,0) ## zero terminate
  }
}
#--------------------------------------------------------------------------------

getxy <- function(dettype, capthist) {
    xy <- xy(capthist)
    ## start[z] indexes the first row in xy 
    ## for each possible count z (including zeros), where z is w-order (isk) 
    start <- abs(capthist)
    start <- head(cumsum(c(0,start)),length(start))
    list(xy = xy, start = start)
}
#--------------------------------------------------------------------------------

recodebinomN <- function (dettype, binomN) {
  binomN <- expandbinomN(binomN, dettype)
  detectr <- names(dettype)  # seems to work 2024-01-29
  detectr[(detectr %in% c('polygon', 'transect')) & (binomN == 0)] <- "poissoncount"
  detectr[(detectr %in% c('polygon', 'transect')) & (binomN > 0)]  <- "binomialcount"
  newbinomN <- function (det,N) {
      switch(det, polygonX = -2, transectX = -2, poissoncount = 0, binomialcount = N, -9)
  }
  out <- mapply(newbinomN, detectr, binomN)
  if (any(out < -9))
    stop("secrpoly not ready for detector type")
  out
}
#--------------------------------------------------------------------------------

nullCH <- function (dimCH, individual) {
    if (is.null(individual)) {
        individual <- TRUE   ## 2020-05-16 for backward compatibility
    }
    if (!individual) {
        dimCH[1] <- 1
    }
    array(0, dim = dimCH)
}
#--------------------------------------------------------------------------------

# compress ch (n x k) by animal: list of detectors with non-zero counts, terminated by -1
# (n x k x 2)
nk2 <- function(ch) {
    one <- function(i) {
        wh <- which(ch[i,1,]>0)
        out <- matrix(-1, nrow = nk, ncol = 2)
        if (length(wh)>0) {
            out[1:length(wh),] <- cbind(ch[i,1,wh], wh-1)
        }
        out
    }
    nk <- dim(ch)[3]
    ## 2022-01-04 catch zero rows
    if (nrow(ch) < 1) {
      array(dim=c(0,nk,2))
    }
    else {
      ch2 <- sapply(1:nrow(ch), one)
      ch3 <- array(ch2, dim = c(nk,2,nrow(ch)))
      aperm(ch3, c(3,1,2))
    }
}
#--------------------------------------------------------------------------------

## 2-D n x s if exclusive detector
## 2-D s x k if capped detector (deferred)
## 3-D n x k x 2 (which detector, count of positive records)
compressCH <- function (CH, binomN, fastproximity) {  
    if (all(binomN == -2)) {
        lost <- apply(CH, 1:2, min)<0
        CH <- abs(CH)
        CH <- apply(CH, 1:2, which.max) *(apply(CH, 1:2, max)>0)
        CH[lost] <- -CH[lost]
    }
    # else if (all(binomN == -3)) {
    #   lost <- apply(CH, 2:3, min)<0
    #   CH <- abs(CH)
    #   CH <- apply(CH, 2:3, which.max) *(apply(CH, 2:3, max)>0)
    #   CH[lost] <- -CH[lost]
    # }
    else if (fastproximity) {
        CH <- nk2(CH) 
    }
    CH
}
############################################################################################

decompressCH <- function (CH, fastproximity) {  
    if (fastproximity) {
        out <- array(0, dim=c(nrow(CH), 1, ncol(CH)))
        for (i in 1:nrow(CH)) {
            n <- which(CH[i,,1]>0)
            k <- CH[i,n,1]
            count <- CH[i,n,2]
            out[i,1,k] <- count
        }
        return(out)
    }
    else {
        return (CH)
    }
}
############################################################################################

##############################################################################

prepareSessionData <- function (capthist, mask, maskusage, 
    design, design0, detectfn, groups, fixed, hcov, details, 
    aslist = TRUE, sessnum = 1) {
    ## aslist used internally to determine whether single-session data are wrapped in a list
    if (ms(capthist)) {
        if (!ms(mask)) stop ("expect session-specific mask in prepareSessionData")
        if (is.null(maskusage))
            maskusage <- vector('list', length(capthist))
        mapply(
            prepareSessionData, 
            capthist = capthist, 
            mask = mask, 
            maskusage = maskusage,
            sessnum = 1:length(capthist),
            MoreArgs = list(design, design0, detectfn, groups, fixed, 
                hcov, details, FALSE), 
            SIMPLIFY = FALSE)
    }
    else {
        nc   <- nrow(capthist)
        s    <- ncol(capthist)
        m    <- nrow(mask)
        traps   <- traps(capthist)
        dettype <- detectorcode(traps, MLonly = TRUE, noccasions = s)
        binomNcode <- recodebinomN(dettype, details$binomN)
        ## k-1 because we have zero-terminated these vectors
        k <- getk(traps)
        K <- if (length(k)>1) length(k)-1 else k
        cumk <- cumsum(c(0,k))[1:length(k)]
        
        ## knownclass for hcov mixture models
        knownclass <- getknownclass(capthist, details$nmix, hcov)
        
        ## get static distance matrix
        distmat2 <- getdistmat2(traps, mask, details$userdist, detectfn == 20)

        n.distrib <- switch (tolower(details$distribution), poisson=0, binomial=1, 0)
        xy <- getxy (dettype, capthist)
        usge <- usage(traps)
        if (is.null(usge) || details$ignoreusage) {
            usge <- matrix(1, nrow = K, ncol = s)
        }
        if (is.null(maskusage)) {
            maskusage <- maskboolean(capthist, mask, details$maxdistance)
        }
        else {
            if (!is.matrix(maskusage) || nrow(maskusage) != nrow(capthist) || ncol(maskusage) != nrow(mask))
                stop ('specified maskusage should be n x m matrix of logical values')
            maskusage[] <- as.logical(maskusage)
        }
        
        if (!is.null(details$externalpdot)) {
            if (!(details$externalpdot %in% names(covariates(mask)))) 
                stop ("externalpdot '", details$externalpdot, "' not found in mask covariates")
            externalpdot <- covariates(mask)[, details$externalpdot]
            message("using external pdot")
        }
        else {
            externalpdot <- NULL
        }
        
        ## Groups
        grp  <- group.factor (capthist, groups)
        if (any(is.na(grp))) {
            stop("group is missing for at least one animal")
        }
        ngroup <- max(1,length(group.levels(capthist, groups)))
        CH <- compressCH(capthist, binomNcode, details$fastproximity) 
        
        # 2023-06-09 tentatively remove hcov from condition
        # this leaves some uncertainty: 
        # when is full CH0 (1 row per animal) really needed?
        # why is this an issue for polygonhistoriescpp and not simplehistoriescpp?
        
        # CH0 <- nullCH(dim(CH), packageVersion('secr')<'4.0.0' || design0$individual || ngroup>1 || !is.null(hcov))   ## all-zero CH
        
        CH0 <- nullCH(dim(CH), packageVersion('secr')<'4.0.0' || design0$individual || ngroup>1)   ## all-zero CH
        
        #####################################################################
        logmult <- 0 # omitted for polygon detectors
        #####################################################################
        
        data <- list(
            sessnum = sessnum,    # added 2021-06-22
            CH = CH,
            CH0 = CH0,
            nc = nc,
            s = s,
            k = k,
            K = K,
            cumk = cumk,
            m = m,
            traps = traps,
            dettype = dettype,
            binomNcode = binomNcode,
            usge = usge,
            mask = mask,
            externalpdot = externalpdot,
            distmat2 = distmat2,
            knownclass = knownclass,
            n.distrib = n.distrib,
            xy = xy,
            grp = grp,
            maskusage = maskusage,
            logmult = logmult,
            HPXpoly = FALSE
        )
    if (aslist) list(data=data)
    else data
    #####################################################################
  }
}
