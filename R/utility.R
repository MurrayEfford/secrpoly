################################################################################
## package 'secrpoly'
## utility.R
## 2024-01-29, 2024-09-20
################################################################################

# Global variables in namespace
#
## define a local environment for temporary variables e.g. iter
.localstuff <- new.env()

.localstuff$packageType <- ' pre-release'
##.localstuff$packageType <- ''

.localstuff$validdetectors <- c('single','multi','proximity','count', 
                                'polygonX', 'transectX', 'signal', 'polygon', 'transect', 
                                'capped', 'null','null','null','null', 'telemetry', 'signalnoise')
.localstuff$polydetectors <- c('polygon','transect','polygonX','transectX')
.localstuff$countdetectors <- c('count','polygon','transect','unmarked','telemetry')
.localstuff$iter <- 0   ## counter 1
.localstuff$iter2 <- 0  ## counter 2
.localstuff$detectionfunctions <-
  c(as.character(0:13),
    'hazard halfnormal',
    'hazard hazard rate',
    'hazard exponential',
    'hazard annular normal',
    'hazard cumulative gamma',
    'hazard variable power')

.localstuff$DFN <- c(as.character(0:13), 
                     'HHN', 'HHR', 'HEX', 'HAN', 'HCG', 'HVP')

.localstuff$learnedresponses <- c('b', 'bk', 'B', 'k', 'Bk') 

#-------------------------------------------------------------------------------

detectionfunctionname <- function (fn) {
    .localstuff$detectionfunctions[fn+1]
}

#-------------------------------------------------------------------------------

detectionfunctionnumber <- function (detname) {
    dfn <- match (toupper(detname), .localstuff$DFN)
    if (is.na(dfn))
        dfn <- match (tolower(detname), .localstuff$detectionfunctions)
    if (is.na(dfn))
        stop ("unrecognised detection function ", detname)
    dfn-1
}

#-------------------------------------------------------------------------------

valid.model <- function(model, CL, detectfn, hcov, userdist, sessioncovnames) {
    if (any(sapply(model, badsmooths)))
        warning ("smooth term may be unsuitable for secr: ",
                 "does not specify k or fx where required")
}

#-------------------------------------------------------------------------------

valid.pnames <- function (details, CL, detectfn, alltelem, sighting, nmix) {
    ## modelled parameters
    pnames <- list(
        c('lambda0','sigma'),      # 14 hazard halfnormal
        c('lambda0','sigma','z'),  # 15 hazard hazard rate
        c('lambda0','sigma'),      # 16 hazard exponential
        c('lambda0','sigma','w'),  # 17
        c('lambda0','sigma','z'),  # 18
        c('lambda0','sigma','z')   # 19
        )[[detectfn-13]]

    if (details$param %in% c(2,6))
        pnames[1] <- 'esa'
    if (details$param %in% c(3,5))
        pnames[1] <- 'a0'
    if (details$param %in% 4:6) {
        pnames[2] <- 'sigmak'
        pnames <- c(pnames, 'c')
        pnames <- c(pnames, 'd')
    }
    if (!CL)
      pnames <- c('D', pnames)
    if (nmix>1)
        pnames <- c(pnames, 'pmix')
    pnames
}
#-------------------------------------------------------------------------------

new.param <- function (details, model, CL) {
    esa <- 'esa' %in% names(model)
    a0 <- 'a0' %in% names(model)
    sigmak <- 'sigmak' %in% names(model)
    newparam <- details$param
    if (esa & !sigmak) {
        newparam <- 2
    }
    if (a0 & !sigmak) {
        newparam <- 3
    }
    if (sigmak) {
        if (esa) {
            newparam <- 6
        }
        else {
            if (CL)
                stop ("sigmak parameterization requires full model, not CL, unless also 'esa'")
            newparam <- ifelse(a0, 5, 4)
        }
    }
    if (newparam  != details$param)
        warning ("Using parameterization details$param = ", newparam, call. = FALSE)
    newparam
}

#-------------------------------------------------------------------------------

expanddet <- function(CH) {
    trps <- traps(CH)
    if (is.null(trps))
        return ('nonspatial')
    else {
        det <- detector(trps)
        if (length(det)<ncol(CH))
            rep(det[1], ncol(CH))
        else det
    }
}

#-------------------------------------------------------------------------------

ndetector <- function (traps) {
    if (is.null(traps))
        return(1)
    else if (all(detector(traps) %in% .localstuff$polydetectors))
        length(levels(polyID(traps)))
    else
        nrow(traps)
}

#-------------------------------------------------------------------------------

memo <- function (text, trace) {
    ## could use message(text), but does not immediately flush console
    if (trace) { cat (text, '\n')
    flush.console() }
}

#-------------------------------------------------------------------------------

## miscellaneous functions

sine     <- function (x) asin (x*2-1)
invsine  <- function (y) (sin(y)+1) / 2
odds     <- function (x) x / (1-x)
invodds  <- function (y) y / (1+y)

#-------------------------------------------------------------------------------

lnbinomial <- function (x,size,prob) {
    # dbinom allowing non-integer x, forcing log = TRUE
    if (x <= size) {
        lgamma (size+1) - lgamma (size-x+1) - lgamma (x+1) +
            x * log(prob) + (size-x) * log (1-prob)
    }
    else {
        -Inf
    }
}

#-------------------------------------------------------------------------------

var.in.model <- function(v,m) v %in% unlist(lapply(m, all.vars))

#-------------------------------------------------------------------------------

get.nmix <- function (model, capthist, hcov) {
    model$D <- NULL  ## ignore density model
    model$pmix <- NULL ## pmix alone cannot make this a mixture model
    nmix <- 1
    if (any(var.in.model('h2', model))) {
        nmix <- 2
        if (any(var.in.model('h3', model)))
            stop ("do not combine h2 and h3")
    }
    if (any(var.in.model('h3', model))) {
        nmix <- 3
    }
    if ((nmix == 1) & (!is.null(hcov))) {
        if (ms(capthist))
            capthist <- capthist[[1]]
        if (is.factor(covariates(capthist)[,hcov]))
            lev <- levels(covariates(capthist)[,hcov])
        else
            lev <- levels(factor(covariates(capthist)[,hcov]))
        if (all(is.na(covariates(capthist)[,hcov])))
            stop ("hcov missing for all individuals, but detection model invariant")
        if (length(lev) < 2)
            stop ("hcov covariate not found or has fewer than 2 levels")
        if (length(lev) > 2)
            warning ("hcov covariate has more than 2 levels; using only first two", call. = FALSE)
        nmix <- 2
    }
    nmix
}

#-------------------------------------------------------------------------------

spatialscale <- function (object, detectfn, session = '') {
    if (inherits(object, 'secr')) {
        if (ms(object))
            detpar <- detectpar(object)[[session]]
        else
            detpar <- detectpar(object)
        cutval <- object$details$cutval
    }
    else {
        detpar <- object
        cutval <- object$cutval
    }
    if (!is.null(detpar$sigma)) detpar$sigma
    else if (detectfn == 10) {
        (cutval - detpar$beta0) / detpar$beta1
    }
    else if (detectfn == 11) {
        d11 <- function(d, beta0, beta1, c) beta0 +
            beta1 * (d-1) - 10 * log10(d^2) - c
        interval <- c(0,10 * (cutval - detpar$beta0) / detpar$beta1)
        uniroot (d11, interval, detpar$beta0, detpar$beta1, cutval)$root
    }
    else if (detectfn == 9) {
        - 1 / detpar$b1   
    }
    else stop ("unrecognised detectfn")
}

#-------------------------------------------------------------------------------

## logical for whether object specifies userDfn
userD <- function (object) {
  if (!inherits(object, c('secr','ipsecr')))
    stop ("requires fitted model")
  !is.null(object$details$userDfn)
}

#-------------------------------------------------------------------------------

## Detection functions

HHN <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]
    1 - exp(-lambda0 * exp (-r^2 / 2 / sigma^2))
}
HHR <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    1 - exp(-lambda0 * ( 1 - exp (-(r / sigma)^-z)))
}
HEX <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]
    1 - exp(-lambda0 * exp (-r / sigma))
}
HAN <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
    1 - exp(-lambda0 * exp (-(r-w)^2 / 2 / sigma^2))
}
HCG <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    lambda0 * pgamma(r, shape=z, scale=sigma/z, lower.tail = FALSE)
}
HVP <- function (r, pars, cutval) {
    lambda0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
    1 - exp(-lambda0 * exp(-(r/sigma)^z))
}

#-------------------------------------------------------------------------------


# transformation tidy up 2021-12-16
# arbitrary link function specified with functions X, invX, se.invX

transform <- function (x, link) {
    switch (link,
        identity = x,
        i1000 = x * 1000,
        log = log(x),
        neglog = log(-x),
        logit = logit(x),
        odds = odds(x),
        sin = sine(x),
        do.call(link, list(x))
    )
}
#-------------------------------------------------------------------------------

untransform <- function (beta, link) {
    switch (link,
        identity = beta,
        i1000 = beta / 1000,
        log = exp(beta),
        neglog = -exp(beta),
        logit = invlogit(beta),
        odds = invodds(beta),
        sin = invsine(beta),
        do.call(paste0('inv',link), list(beta))
    )
}
#-------------------------------------------------------------------------------

mlogit.untransform <- function (beta, latentmodel) {
    if (!missing(latentmodel)) {
        for (i in unique(latentmodel))
            beta[latentmodel==i] <- mlogit.untransform(beta[latentmodel==i])
        beta
    }
    else {
        ## beta should include values for all classes (mixture components)
        nmix <- length(beta)
        if (sum(is.na(beta)) != 1) {
            ## require NA for a single reference class
            rep(NA, length(beta))
        }
        else {
            nonreference <- !is.na(beta)   # not reference class
            b <- beta[nonreference]
            pmix <- numeric(nmix)
            pmix[nonreference] <- exp(b) / (1+sum(exp(b)))
            pmix[!nonreference] <- 1 - sum(pmix[nonreference])
            pmix
        }
    }
}

#-------------------------------------------------------------------------------

mlogit <- function (x) {
    ## return the mlogit of an unscaled vector of positive values
    logit(x/sum(x))
}

## End of miscellaneous functions

#-------------------------------------------------------------------------------

group.levels <- function (capthist, groups, sep='.') {
    # 2016-06-05 use also for trap strata
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, group.levels, groups, sep)
        sort(unique(unlist(temp)))  ## vector of global levels
    }
    else {
        if (is.null(groups)) 0
        else {
            if (!all(groups %in% names(covariates(capthist))))
                stop ("one or more grouping variables is missing ",
                      "from covariates")
            temp <- as.data.frame(covariates(capthist)[,groups])
            # omit null combinations, sort as with default of factor levels
            sort(levels(interaction(temp, drop=T, sep=sep)))
        }
    }
}

#-------------------------------------------------------------------------------


group.factor <- function (capthist, groups, sep='.')
    ## convert a set of grouping factors to a single factor (g)
    ## levels common to all sessions
{
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, group.factor, groups)  ## recursive call
        grouplevels <- group.levels(capthist, groups)
        if (length(grouplevels)<2)
            temp
        else
            # list; force shared factor levels on each component
            lapply (temp, factor, levels=grouplevels)
    }
    else {
        if (is.null(groups) | (length(groups)==0) )
            return (factor(rep(1, nrow(capthist)), levels = 1))  # added levels 2017-04-18
        temp <- as.data.frame(covariates(capthist)[,groups])
        if (ncol(temp) != length(groups))
            stop ("one or more grouping variables is missing from ",
                  "covariates(capthist)")
        temp <- interaction(temp, drop=T, sep=sep)  # omit null combinations
        temp
    }
}

#-------------------------------------------------------------------------------

## Return an integer vector of class membership defined by a categorical
## individual covariate in a capthist object. Individuals of unknown
## class (including those with class exceeding nmix) are coded 1,
## others as (class number + 1). When no mixture is specified (nmix == 1)
## all are coded as unknown.

## knownclass 1 'unknown' 
## knownclass 2 'latent class 1' 
## knownclass 3 'latent class 2' 

getknownclass <- function(capthist, nmix, hcov) {
    if (ms(capthist)) {
        lapply(capthist, getknownclass, nmix = nmix, hcov = hcov)
    }
    else {
        if ((nmix>1) & (!is.null(hcov))) {
          ## 2020-09-05 use as.factor() instead of factor() to coerce 
          ## (if already factor, coercing with factor() loses old levels)
          var <- as.factor(covariates(capthist)[,hcov])
          tmp <- as.numeric(var) + 1
          tmp[is.na(tmp) | (tmp>(nmix+1))] <- 1
          attr(tmp,'levels') <- levels(factor(covariates(capthist)
            [,hcov]))[1:nmix]
          tmp
        }
        else
            rep(1,nrow(capthist))
    }
}

#-------------------------------------------------------------------------------

getnmix <- function (details) {
    if (is.null(details$nmix))
       1
    else
       details$nmix
}

#-------------------------------------------------------------------------------

## expand beta parameter vector using template of 'fixed beta'
## fixed beta fb input is missing (NA) for estimated beta parameters
fullbeta <- function (beta, fb) {
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta  ## partial beta (varying only)
        beta <- fb             ## complete beta
    }
    beta
}

#-------------------------------------------------------------------------------

getbinomN <- function (binomN, detectr) {
    if (any(detectr %in% .localstuff$countdetectors)) {
        if (is.null(binomN))
            return(0)
        else if (any(binomN == 'usage'))
            return(1)
        else
            return(binomN)
    }
    else
        return(1)
}

#-------------------------------------------------------------------------------

complete.beta <- function (object) {
    fb <- object$details$fixedbeta
    # modified 2022-04-02 for consistency with ipsecr
    beta <- if (inherits(object, 'secr')) object$fit$par else object$beta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        fb[is.na(fb)] <- beta
        beta <- fb
    }
    beta
}

#-------------------------------------------------------------------------------

complete.beta.vcv <- function (object) {
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        beta.vcv <- matrix(NA, nrow = nbeta, ncol = nbeta)
        beta.vcv[is.na(fb[row(beta.vcv)]) & is.na(fb[col(beta.vcv)])] <- object$beta.vcv
    }
    else {
        beta.vcv <- object$beta.vcv
    }
    beta.vcv
}

#-------------------------------------------------------------------------------

smooths <- function (formula) {
    ## which terms in formula are smooths?
    ## returns logical vector
    labels <- attr(terms(formula), "term.labels")
    if (length(labels) > 0)
        sapply(labels, function (x) any(sapply(c("s\\(", "te\\(", "poly\\("), grepl, x)))
    else
        logical(0)
}

#-------------------------------------------------------------------------------

polys <- function (formula) {
    ## which terms in formula are orthogonal polynomials?
    ## returns logical vector
    labels <- attr(terms(formula), "term.labels")
    if (length(labels) > 0)
        sapply(labels, grepl, pattern = "poly\\(")
    else
        logical(0)
}

#-------------------------------------------------------------------------------

badsmooths <- function (formula) {
    ## does smooth specification conform to secr requirements?
    ## returns TRUE/FALSE
    labels <- attr(terms(formula), "term.labels")
    if (length(labels) > 0) {
        smoothterms <- sapply(labels, function (x)
                              any(sapply(c("s\\(", "te\\("), grepl, x)))
        labels <- labels[smoothterms]
        any(sapply(labels, function(x)
               grepl("s\\(", x) & !grepl("k =", x))) |
        any(sapply(labels, function(x)
               grepl("te\\(", x) & (!grepl("fx = TRUE", x) | !grepl("k =", x))))
    }
    else
        FALSE
}

#-------------------------------------------------------------------------------

gamsetup <- function(formula, data, ...) {
    ## use 'session' column as dummy LHS so gam does not gag
    ## (cf secrgam:::make.density.design.matrix)
    ## session is always present in detection data, must be added for D
    if (is.null(data$session)) data$session <- rep(1,nrow(data))
    formula <- update.formula(formula, session ~ .)
    setup <- gam(formula, data = data, fit = FALSE, ...)
    colnames(setup$X) <- setup$term.names
    setup
}
#-------------------------------------------------------------------------------

makerealparameters <- function (design, beta, parindx, link, fixed) {
    modelfn <- function(i) {
        ## linear predictor for real parameter i
        Yp <- design$designMatrices[[i]] %*% beta[parindx[[i]]]
        if (names(link)[i] == 'pmix') {
            ## 2013-04-14 index of class groups (pmix sum to 1.0 within latentmodel)
            cols <- dimnames(design$designMatrices[[i]])[[2]]
            h2 <- grep('.h2', cols, fixed=T)
            h3 <- grep('.h3', cols, fixed=T)
            h2c <- grep(':h2', cols, fixed=T)
            h3c <- grep(':h3', cols, fixed=T)
            h.cols <- c(h2,h3,h2c,h3c)
            tmp <- design$designMatrices[[i]][,-h.cols, drop = FALSE]
            tmph <- design$designMatrices[[i]][,h.cols, drop = FALSE]
            ## 2018-02-23 why as.numeric()? 
            latentmodel <- as.numeric(factor(apply(tmp,1,paste, collapse='')))
            refclass <- apply(tmph,1,sum) == 0
            Yp[refclass] <- NA
            Yp <- mlogit.untransform(Yp, latentmodel)
            Yp[design$parameterTable[,i]]
        }
        else {
            Yp <- untransform(Yp, link[[i]])
            Yp[design$parameterTable[,i]]   ## replicate as required
        }
    }
    ## construct matrix of detection parameters
    nrealpar  <- length(design$designMatrices)
    parindx$D <- NULL ## detection parameters only
    link$D    <- NULL ## detection parameters only
    parindx$noneuc <- NULL ## detection parameters only
    link$noneuc    <- NULL ## detection parameters only
    detectionparameters <- names(link)
    fixed.dp <- fixed[detectionparameters[detectionparameters %in% names(fixed)]]
    
    if (length(fixed.dp)>0)
        for (a in names(fixed.dp))  ## bug fixed by adding this line 2011-09-28
            link[[a]] <- NULL
    if (length(link) != nrealpar)
        stop ("number of links does not match design matrices")
    
    if (nrealpar == 0) {
        return(matrix(unlist(fixed.dp),nrow = 1))
    }
    
    temp <- sapply (1:nrealpar, modelfn)
    if (nrow(design$parameterTable)==1) temp <- t(temp)
    nrw <- nrow(temp)
    ## make new matrix and insert columns in right place
    temp2 <- as.data.frame(matrix(nrow = nrw, ncol = length(detectionparameters)))
    names(temp2) <- detectionparameters
    temp2[ , names(design$designMatrices)] <- temp          ## modelled
    if (!is.null(fixed.dp) & length(fixed.dp)>0)
        temp2[ , names(fixed.dp)] <- sapply(fixed.dp, rep, nrw)    ## fixed
    as.matrix(temp2)
    
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

getcellsize <- function (mask) {
    if (inherits(mask, 'linearmask'))
        cell <- attr(mask, 'spacing') / 1000  ## per km
    else
        cell <- attr(mask, 'area')            ## per ha
    if (is.null(cell))
        stop ("mask lacks valid cell size (area or spacing)")
    cell
}

#-------------------------------------------------------------------------------

expandbinomN <- function (binomN, detectorcodes) {
    # assumes detectorcodes is a vector of length = noccasions
    binomN <- ifelse (detectorcodes %in% c(2,6,7), binomN, 1)
    if (any(is.na(binomN))) stop ("NA value in binomN")
    binomN
}

#-------------------------------------------------------------------------------

allzero <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, allzero)
    }
    else {
        telemocc <- detector(traps(object))=='telemetry'
        apply(object[,!telemocc,,drop=FALSE],1,sum)==0
    }
}

#-------------------------------------------------------------------------------

## return indices of first occasion and detector for which PIAx is non-zero 
firstsk <- function (PIAx) {
  ## PIAx dim n,s,k
  wh <- function(d2) {
    match(TRUE, d2>0)
  }
  apply(PIAx,1,wh)
}

#-------------------------------------------------------------------------------

selectCHsession <- function(capthist, sessnum) {
    if (ms(capthist)) 
        capthist[[sessnum]]
    else 
        capthist
}

#-------------------------------------------------------------------------------

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

telemcode <- function(object, ...) {
    if (inherits(object, 'traps') && !ms(object))
        switch (telemetrytype(object), none = 0, 
                independent = 1, dependent = 2, concurrent = 3, 0)
    else 
        NA
}

#-------------------------------------------------------------------------------

recodebinomN <- function (dettype, binomN, telemcode) {
    binomN <- expandbinomN(binomN, dettype)
    detectr <- .localstuff$validdetectors[dettype+2]
    detectr[(detectr %in% c('count','polygon', 'transect')) & (binomN == 0)] <- "poissoncount"
    detectr[(detectr %in% c('count','polygon', 'transect')) & (binomN > 0)]  <- "binomialcount"
    newbinomN <- function (det,N) {
        recoded <- switch(det, 
                          single = -2, multi = -2, polygonX = -2, transectX = -2,
                          proximity = -1, signal = -1, capped = -1, 
                          poissoncount = 0, binomialcount = N, telemetry = -3, -9)
        ## dummy value -9 indicates no action
        if (recoded == -3 && telemcode == 0) recoded <- -7
        recoded
    }
    out <- mapply(newbinomN, detectr, binomN)
    if (any(out < -9))
        stop("secrpoly not ready for detector type")
    out
}
#--------------------------------------------------------------------------------

getpID <- function(PIA, realparval, MRdata)
{
    ss <- dim(PIA)[3]
    nmix <- dim(PIA)[5]
    pID <- matrix(1, nrow = ss, ncol = nmix)
    if ('pID' %in% colnames(realparval)) {
        nc <- dim(PIA)[2]
        if (!is.null(MRdata$Tm) || !is.null(MRdata$Tn)) {
            for (s in 1:ss) {
                if (MRdata$markocc[s]<1)
                    for (x in 1:nmix) {
                        k <- match(TRUE, PIA[1,1,s,,x]>0)
                        c <- PIA[1,1,s,k,x]
                        pID[s,x] <- realparval[c, 'pID'] 
                    }
            }
        }
    }
    pID
}
#--------------------------------------------------------------------------------

## mixture proportions by animal        
## assume dim(PIA)[1] == 1
getpmix <- function(knownclass, PIA, realparval)
{
    nc <- dim(PIA)[2]
    # not needed nc <- length(knownclass)   ## 2020-11-04
    k <- dim(PIA)[4]
    nmix <- dim(PIA)[5]
    pmixn <- matrix(1, nrow = nmix, ncol = nc)
    pmix <- numeric(nmix)
    if (nmix>1) {
        # index of first non-missing occasion s and detector k
        fsk <- sapply(1:nc, function(i) firstsk(PIA[1,i,,,1, drop = FALSE]))
        kc <- as.vector((fsk-1) %/% k + 1)
        sc <- as.vector((fsk-1) %/% k + 1)
        for (x in 1:nmix) {
            c <- PIA[cbind(1,1:nc,sc,kc,x)]
            pmixx <- realparval[c, 'pmix']    ## NOT CONSISTENT WITH pmix numeric(nmix)
            ## knownclass=2 maps to x=1 
            pmixn[x,] <- ifelse (knownclass > 1,
                                 ifelse (knownclass == (x+1), 1, 0),
                                 pmixx)
        }
        ## need pmix for each group... not ready yet
        attr(pmixn, 'pmix') <-  realparval[PIA[cbind(1,1,sc[1],kc[1],1:nmix)],'pmix']
    }
    pmixn
}
#--------------------------------------------------------------------------------

# previously in preparedata.R

getxy <- function(dettype, capthist) {
    xy <- xy(capthist)
    ## start[z] indexes the first row in xy 
    ## for each possible count z (including zeros), where z is w-order (isk) 
    start <- abs(capthist)
    start <- head(cumsum(c(0,start)),length(start))
    list(xy = xy, start = start)
}
#--------------------------------------------------------------------------------
