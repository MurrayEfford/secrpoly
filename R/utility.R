################################################################################
## package 'secrpoly'
## utility.R
## 2024-01-29
################################################################################

# Global variables in namespace
#
## define a local environment for temporary variables e.g. iter
.localstuff <- new.env()

.localstuff$packageType <- ' pre-release'
##.localstuff$packageType <- ''

.localstuff$validdetectors <- c('polygonX', 'transectX', 'polygon', 'transect')
.localstuff$individualdetectors <- c('polygonX', 'transectX', 'polygon', 'transect')
.localstuff$polydetectors <- c('polygon','transect','polygonX','transectX')
.localstuff$exclusivedetectors <- c('polygonX','transectX')
.localstuff$countdetectors <- c('polygon','transect')
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
                     'HHN', 'HHR', 'HEX', 'HAN', 'HCG', 'HVP','HPX')

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

parnames <- function (detectfn) {
    list(
        c('lambda0','sigma'),
        c('lambda0','sigma','z'),
        c('lambda0','sigma'),
        c('lambda0','sigma','w'),
        c('lambda0','sigma','z'),
        c('lambda0','sigma','z')
        )[[detectfn-13]]
}

#-------------------------------------------------------------------------------

valid.detectfn <- function (detectfn, valid = 14:19) {
    if (is.null(detectfn))
        stop ("requires 'detectfn'")
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if (!(detectfn %in% valid))
        stop ("invalid detection function")
    detectfn
}

#-------------------------------------------------------------------------------

valid.model <- function(model, CL, detectfn, hcov, userdist, sessioncovnames) {
    if (any(sapply(model, badsmooths)))
        warning ("smooth term may be unsuitable for secr: ",
                 "does not specify k or fx where required")
}

#-------------------------------------------------------------------------------

getuserdistnames <- function (userdist) {
    ## return the names of any supplementary arguments of user-provided function
    ## for non-euclidean distance computations
    if (is.function(userdist)) {
        distnames <- try(userdist(), silent = TRUE)
        if (!is.character(distnames))
            stop("invalid userdist function - ",
                 "should return parameter names when called with no arguments")
        distnames
    }
    else
        character(0)
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

valid.userdist <- function (userdist, detector, xy1, xy2, mask, sessnum) {
    if (!is.null(userdist)) 
        stop ("userdist cannot be used with polygon detector types;")
    edist(xy1, xy2)      ## Euclidean distance
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

detectorcode <- function (object, MLonly = TRUE, noccasions = NULL) {
    ## numeric detector code from a traps object
    detcode <- sapply(detector(object), switch,
        polygonX    = 3,
        transectX   = 4,
        polygon     = 6,
        transect    = 7,
        -2)
    
    if (MLonly) {
        detcode <- ifelse (detcode==-1, rep(0,length(detcode)), detcode)
        if (any(detcode<0))
            stop ("Unrecognised detector type")
    }

    if (!is.null(noccasions) & (length(detcode)==1))
        detcode <- rep(detcode, noccasions)
    detcode
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

ndetectpar <- function (detectfn) {
    length(parnames(detectfn))
}

#-------------------------------------------------------------------------------

replacedefaults <- function (default, user) replace(default, names(user), user)

#-------------------------------------------------------------------------------

discreteN <- function (n, N) {
    tN <- trunc(N)
    if (N != tN) tN + sample (x = c(1,0), prob = c(N-tN, 1-(N-tN)),
        replace = T, size = n)
    else rep(tN,n)
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

pad1 <- function (x, n) {
## pad x to length n with dummy (first value)
    if (is.factor(x)) {
        xc <- as.character(x)
        xNA <- c(xc, rep(xc[1], n-length(xc)))
        out <- factor(xNA, levels=levels(x))
    }
    else out <- c(x, rep(x[1], n-length(x)))
    out
}

#-------------------------------------------------------------------------------

padarray <- function (x, dims) {
    temp <- array(dim=dims)
    dimx <- dim(x)
    if (all(dimx>0)) {
        if (length(dimx)<2 | length(dimx)>3)
            stop ("invalid array")
        if (length(dimx)>2) temp[1:dimx[1], 1:dimx[2], 1:dimx[3]] <- x
        else temp[1:dimx[1], 1:dimx[2]] <- x
    }
    temp
}

#-------------------------------------------------------------------------------

## regularize a list of formulae
stdform <- function (flist) {
    LHS <- function (form) {
        trms <- as.character (form)
        if (length(trms)==2) '' else trms[2]
    }
    RHS <- function (form) {
        trms <- as.character (form)
        ## 2020-05-14 for compatibility with R 4.0
        if (length(trms)==3) as.formula(paste(trms[c(1,3)], collapse = " ")) else form
    }
    lhs <- sapply(flist, LHS)
    temp <- lapply(flist, RHS)
    if (is.null(names(flist))) names(temp) <- lhs
    else names(temp) <- ifelse(names(flist) == '', lhs, names(flist))
    temp
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

model.string <- function (model, userDfn) {
    # 2023-04-16 Note: model should be a list
    if (!is.null(userDfn)) {
        if (!is.null(model$D))
            model$D <- paste('~userD', userDfn('name'), sep='.')
    }
    temp <- paste (names(model), as.character(model), collapse=' ', sep='')
    temp
}

#-------------------------------------------------------------------------------

fixed.string <- function (fixed) {
    if (is.null(fixed) | length(fixed)==0) 'none'
    else paste (names(fixed), as.character(fixed), collapse=', ', sep=' = ')
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

fixpmix <- function(x, nmix) {

    ## x is a list with component pmix that is a matrix (dataframe)
    ## with columns 'estimate' and 'se' (and possibly others)
    ## referring to the linear predictor of pmix (i.e. on mlogit
    ## scale) and rows corresponding to rows in newdata
    ## (i.e. arbitrary combinations of predictors, including mixture
    ## class h2 or h3)

    ####################################################
    ## It is necessary that newdata include all levels
    ## of the mixture class.
    ####################################################

    ## 2013-10-29
    ## assuming mixture is always last dimension...

    ## previously used in collate, model.average and predict.secr
    ## 2015-09-30 incorporated in secr.lpredictor

    temp <- matrix(x$pmix[,'estimate'], ncol = nmix)
    if (nmix==2) temp[,x$pmix[,'h2']] <- x$pmix[,'estimate']
    if (nmix==3) temp[,x$pmix[,'h3']] <- x$pmix[,'estimate']
    temp2 <- apply(temp, 1, clean.mlogit)
    x$pmix[,'estimate'] <- as.numeric(t(temp2))
    if (nmix==2)
        x$pmix[as.numeric(x$pmix$h2)==1,'se'] <- x$pmix[as.numeric(x$pmix$h2)==2,'se']
    else
        x$pmix[,'se'] <- rep(NA, nrow(x$pmix))   ## don't know how
    x
}

#-------------------------------------------------------------------------------

add.cl <- function (df, alpha, loginterval, lowerbound = 0) {

## add lognormal or standard Wald intervals to dataframe with columns
## 'estimate' and 'SE.estimate'
## lowerbound added 2011-07-15
    z <- abs(qnorm(1-alpha/2))
    if (loginterval) {
        delta <- df$estimate - lowerbound
        df$lcl <- delta / exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
        df$ucl <- delta * exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
    }
    else {
        df$lcl <- pmax(lowerbound, df$estimate - z * df$SE.estimate)
        df$ucl <- df$estimate + z * df$SE.estimate
    }
    df
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


nclusters <- function (capthist) {
    if (ms(capthist)) {
	lapply(capthist, nclusters)
    }
    else 	{
        nmash <- attr(capthist, 'n.mash')
        ifelse (is.null(nmash), 1, length(nmash))
    }
}

#-------------------------------------------------------------------------------

## clunky but effective re-write 2012-09-04, improved 2016-02-20, 2016-05-10
leadingzero <- function (x) {
    xc <- as.character(x)
    w <- max(nchar(xc))
    n0 <- function(n) paste(rep('0',n), collapse='')
    paste(sapply(w-nchar(xc), n0), x, sep='')

    ## or, 2016-01-15, 2016-02-20 BUT DOESN'T HANDLE NON-INTEGER 2016-05-10
    #     if (is.character(x)) x <- as.numeric(x)
    #     sprintf(paste("%0", w, "d", sep = ""), x)
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

gradient <- function (pars, fun, eps=0.001, ...)
## quick & dirty 2009 09 14
## used by plot.secr for delta method limits
{
  est <- pars
  g   <- pars
  for (i in 1:length(est))
  {
      temp     <- est[i]
      if (temp != 0.0) delta <- eps * abs(temp)
      else             delta <- eps
      est[i]  <- temp - delta
      fminus  <- fun (est, ...)
      est[i]  <- temp + delta
      fplus   <- fun (est, ...)
      g[i]    <- (fplus - fminus) / (2.0 * delta)
      est[i]  <- temp;
  }
  g
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

# used only in model.average, modelAverage
se.transform <- function (real, sereal, link) {
    switch (link,
        identity = sereal,
        i1000 = sereal / 1000,
        log = log((sereal/real)^2 + 1)^0.5,
        neglog = log((sereal/-real)^2 + 1)^0.5,
        logit = sereal / real / (1 - real),
        sin = NA,
        do.call(paste0('se.',link), list(real, sereal) )
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

se.untransform <- function (beta, sebeta, link) {
    # Approximate translation of SE to untransformed scale
    # Delta method cf Lebreton et al 1992 p 77
    switch (link,
        identity = sebeta,
        i1000 = sebeta / 1000,
        log = exp(beta) * sqrt(exp(sebeta^2)-1),
        neglog = exp(beta) * sqrt(exp(sebeta^2)-1),
        logit = invlogit(beta) * (1-invlogit(beta)) * sebeta,
        sin = NA,                ####!!!!
        do.call(paste0('se.inv', link), list(beta=beta, sebeta=sebeta))
    )
}
#-------------------------------------------------------------------------------

# vectorized transformations
Xtransform <- function (real, linkfn, varnames) {
    mapply(transform, real, linkfn[varnames])
}
se.Xtransform <- function (real, sereal, linkfn, varnames) {
    mapply(se.transform, real, sereal, linkfn[varnames])
}
Xuntransform <- function (beta, linkfn, varnames) {
    mapply(untransform, beta, linkfn[varnames])
}
se.Xuntransform <- function (beta, sebeta, linkfn, varnames)
{
    if (length(beta)!=length(sebeta))
        stop ("'beta' and 'sebeta' do not match")
    if (!all(varnames %in% names(linkfn)))
        stop ("'linkfn' component missing for at least one real variable")
    mapply(se.untransform, beta, sebeta, linkfn[varnames])
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

clean.mlogit <- function(x) {
    ## 2014-08-19 for robustness...
    if (is.na(x[2])) x[2] <- 1-x[1]
    x[1] <- NA   ## assumed reference class
    logit(mlogit.untransform(x))
}

#-------------------------------------------------------------------------------

mlogit <- function (x) {
    ## return the mlogit of an unscaled vector of positive values
    ## 2013-04-14
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

h.levels <- function (capthist, hcov, nmix) {
    ## determine the first nmix levels of a factor individual covariate
    if (is.null(hcov))
        as.character(1:nmix)
    else {
        if (ms(capthist)) {
            ## take first session as we can assume factor covariates have same levels in
            ## all sessions
            capthist <- capthist[[1]]
        }
        hcov <- covariates(capthist)[,hcov]
        if (!is.factor(hcov)) {
            warning ("hcov was coerced to a factor", call. = FALSE)
            hcov <- factor(hcov)
        }
        levels(hcov)[1:nmix]
    }
}

#-------------------------------------------------------------------------------

n.occasion <- function (capthist) {
## return the number of sampling occasions for each session in capthist
    if (inherits(capthist, 'list')) {
        sapply(capthist, n.occasion)
    }
    else {
        ncol(capthist)
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

getgrpnum <- function (capthist, groups) {
    if (is.null(groups))
        rep(1, nrow(capthist))
    else
        match(group.factor(capthist, groups), group.levels(capthist, groups))
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

## inflate a convex outline along all radii by linear factor 'rmult'
inflate <- function (xy, rmult = 1) {
    xy <- as.matrix(xy)
    centre <- apply(xy, 2, mean)
    xy <- sweep(xy, MARGIN = 2, STATS = centre, FUN = '-')
    r <- apply(xy, 1, function(z) sqrt(sum(z^2)))
    theta <- atan2 (xy[,2], xy[,1])
    r <- r * rmult
    xy <- cbind(r * cos(theta), r * sin(theta))
    sweep(xy, MARGIN = 2, STATS = centre, FUN = '+')
}

#-------------------------------------------------------------------------------

## moved from pdot.R 2013-11-09
## scalar 2016-10-14
getbinomN <- function (binomN, detectr) {
    if (any(detectr %in% .localstuff$countdetectors)) {
        if (is.null(binomN))
            return(0)
        else if (binomN == 'usage')
            return(1)
        else
            return(binomN)
    }
    else
        return(1)
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

inflatechull <- function (poly, r, ntheta = 60) {
    theta <- (2*pi) * (1:ntheta) / ntheta
    ## add supernumerary vertices
    temp  <- data.frame(x = apply(expand.grid(poly$x, r * cos(theta)),1,sum),
                   y = apply(expand.grid(poly$y, r * sin(theta)),1,sum))
    hull <- chull(temp)
    temp[c(hull,hull[1]), ]
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

secr.lpredictor <- function (formula, newdata, indx, beta, field, beta.vcv=NULL,
    smoothsetup = NULL, contrasts = NULL, Dfn = NULL) {
    ## form linear predictor for a single 'real' parameter
    ## smoothsetup should be provided whenever newdata differs from
    ## data used to fit model and the model includes smooths from gam
    vars <- all.vars(formula)
    OK <- vars %in% names(newdata)
    if (any(!OK)) {
        missingvars <- paste(vars[!OK], collapse = ', ')
        if (sum(!OK) == 1)
            stop ("model covariate ", missingvars, " not found in 'newdata'")
        else
            stop ("model covariates ", missingvars, " not found in 'newdata'")
    }
    newdata <- as.data.frame(newdata)
    lpred <- matrix(ncol = 2, nrow = nrow(newdata), dimnames = list(NULL,c('estimate','se')))

    if (!is.null(Dfn) && field == 'D') {
        warning("secr.lpredictor is not ready for D as function -  do not use estimates")
       nsess <- length(unique(newdata$session))
       Yp <- Dfn(newdata[,vars[1]], beta = beta[indx], dimD = c(nrow(newdata)/nsess,1,nsess)) 
       mat <- as.matrix(newdata[,vars[1], drop = FALSE])
    }
    else {
        
        mat <- secr:::general.model.matrix(formula, data = newdata, gamsmth = smoothsetup, 
            contrasts = contrasts)
        if (nrow(mat) < nrow(newdata))
            warning ("missing values in predictors?", call. = FALSE)
        
        nmix <- 1
        if (field=='pmix') {
            ## drop pmix beta0 column from design matrix (always zero)
            mat <- mat[,-1,drop=FALSE]
            if ('h2' %in% names(newdata)) nmix <- 2
            if ('h3' %in% names(newdata)) nmix <- 3
            mixfield <- c('h2','h3')[nmix-1]
        }
        
        ###############################
        Yp <- mat %*% beta[indx]
        ###############################
        
        ## A latent model comprises one row for each latent class.
        ## Back transformation of pmix in mlogit.untransform() requires all rows of 
        ## each latent model. That function splits vector Yp by latent model.
        
        if (field == 'pmix') {
            nonh <- newdata[, names(newdata) != mixfield, drop = FALSE]
            latentmodel <- factor(apply(nonh, 1, paste, collapse = ''))
            refclass <- as.numeric(newdata[, mixfield]) == 1
            Yp[refclass] <- NA
            Yp <- mlogit.untransform(Yp, latentmodel)
            Yp <- logit(Yp)  # return to logit scale for later untransform!
            if (nmix==2) {
                h2.1 <- as.numeric(newdata$h2)==1
                h2.2 <- as.numeric(newdata$h2)==2
            }
        }
    }

    lpred[,1] <- Yp
    if (is.null(beta.vcv) || (any(is.na(beta[indx])))) return ( cbind(newdata,lpred) )
    else {
        if (is.null(Dfn) || field != 'D') {
            vcv <- beta.vcv[indx,indx, drop = FALSE]
            vcv[is.na(vcv)] <- 0
            nrw <- nrow(mat)
            vcv <- apply(expand.grid(1:nrw, 1:nrw), 1, function(ij)
                mat[ij[1],, drop=F] %*% vcv %*% t(mat[ij[2],, drop=F])) 
            
            vcv <- matrix (vcv, nrow = nrw)
            if (field=='pmix') {
                if (nmix==2)
                    vcv[h2.1,h2.1] <- vcv[h2.2,h2.2]
                else
                    vcv[,] <- NA
            }
            lpred[,2] <- diag(vcv)^0.5
        }
        else {
            vcv <- NULL
        }
        
        temp <- cbind(newdata,lpred)
        attr(temp, 'vcv') <- vcv
        return(temp)
    }
}

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

## intercept and fix certain models with bad defaults
updatemodel <- function (model, detectfn, detectfns, oldvar, newvar, warn = FALSE) {
    if (detectfn %in% detectfns) {
        for (i in 1:length(oldvar)) {
            if (oldvar[i] %in% names(model)) {
                names(model)[names(model) == oldvar[i]] <- newvar[i]
                if (warn)
                    warning ("replacing ", oldvar[i], " by ", newvar[i],
                             " in model for detectfn ", detectfn)
            }
        }
    }
    model
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

nparameters <- function (object) {
    Npar <- max(unlist(object$parindx))
    Npar <- Npar + length(object$details$miscparm)
    ## allow for fixed beta parameters
    if (!is.null(object$details$fixedbeta))
        Npar <- Npar - sum(!is.na(object$details$fixedbeta))
    Npar
}

#-------------------------------------------------------------------------------

mapbeta <- function (parindx0, parindx1, beta0, betaindex)

    ## Extend beta vector from simple model (beta0) to a more complex (i.e. general)
    ## model, inserting neutral values (zero) as required.
    ## For each real parameter, a 1:1 match is assumed between
    ## beta values until all beta values from the simpler model are
    ## used up. THIS ASSUMPTION MAY NOT BE JUSTIFIED.
    ## betaindex is a user-controlled alternative.

{
    ## list of zeroed vectors, one per real parameter
    beta1 <- lapply(parindx1, function (x) {x[]<-0; x})

    if (!is.null(betaindex)) {
        beta1 <- unlist(beta1)
        if (sum(betaindex>0) != length(beta0))
            stop ("invalid 'betaindex'")
        beta1[betaindex] <- beta0
        beta1
    }
    else {
        ## indx is within-parameter rather than absolute index
        ## for each _original_ real parameter
        indx <- lapply(parindx0, function(x) x-x[1]+1)
        ## for (j in 1:length(beta1))
        ## improved replace by name2015-11-17
        for (j in names(beta1)) {
            if (j %in% names(beta0))
                beta1[[j]][indx[[j]]] <- beta0[parindx0[[j]]]
        }
        unlist(beta1)
    }
}

#-------------------------------------------------------------------------------

xyinpoly <- function (xy, trps) {
    ptinside <- function (i,k) {
        ## is point i inside poly k?
        polyxy <- as.matrix(lxy[[k]])
        polyxy <- rbind(polyxy, polyxy[1,])   ## close 2014-08-28
        nr <- nrow(polyxy)
        temp <- insidecpp(unlist(xy[i,]), 0, nr-1, as.matrix(polyxy))
    }
    lxy <- split (trps, polyID(trps))
    firstinside <- function (i) {
        frstk <- 0
        for (k in 1:length(lxy)) {
            if (ptinside(i,k)) {
                frstk <- k
                break
            }
        }
        frstk
    }
    sapply(1:nrow(xy), firstinside)
}

#-------------------------------------------------------------------------------

addzerodf <- function (df, oldCH, sess) {
    ## add dummy detection records to dataframe for 'all-zero' case
    ## that arises in sighting-only mark-resight with known marks
    allzero <- apply(oldCH,1,sum)==0
    naz <- sum(allzero)
    if (naz > 0) {
        df0 <- expand.grid(
          newID = rownames(oldCH)[allzero], 
          newocc = NA,
          newtrap = trap(oldCH)[1], 
          alive = TRUE, 
          sess = sess,
          stringsAsFactors = FALSE)
        df$x <- NULL; df$y <- NULL  ## 2021-04-08
        df <- rbind(df,df0)
        if (!is.null(xy(oldCH))) {
            df$x <- c(xy(oldCH)$x, rep(NA, naz))
            df$y <- c(xy(oldCH)$y, rep(NA, naz))
        }
        if (!is.null(signal(oldCH)))  {
            df$signal <- c(signal(oldCH), rep(NA, naz))
        }
    }
    df
}

#-------------------------------------------------------------------------------

expandbinomN <- function (binomN, detectorcodes) {
    # assumes detectorcodes is a vector of length = noccasions
    binomN <- ifelse (detectorcodes %in% c(2,6,7), binomN, 1)
    if (any(is.na(binomN))) stop ("NA value in binomN")
    binomN
}

#-------------------------------------------------------------------------------

newstr <-function (strings) {
    ## compress a character vector
    ## use run length encoding function
    rl <- rle(strings)
    st <- rl$values
    le <- paste0(' (',as.character(rl$lengths), ')')
    le[le==' (1)'] <- ''
    paste(paste0(st, le), collapse = ', ')
}
# newstr(c("single", rep("proximity",4)))

#-------------------------------------------------------------------------------

outsidemask <- function(CH, mask, threshold = spacing(mask) / sqrt(2)) {
    xylist <- telemetryxy(CH)
    dfun <- function(xy) {
        centres <- matrix(apply(xy, 2, mean), ncol = 2)
        distancetotrap(centres, mask)
    }
    sapply(xylist, dfun) > threshold
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

primarysessions <- function(intervals) {
    primarysession <- cumsum(c(0,intervals))
    match(primarysession, unique(primarysession))
}

#-------------------------------------------------------------------------------

secondarysessions <- function(intervals) {
    primary <- primarysessions(intervals)
    unname(unlist(sapply(table(primary), seq_len)))  
}

#-------------------------------------------------------------------------------

boundarytoSF <- function (poly) {
  if (is.null(poly)) {
    NULL
  }
  else if(inherits(poly, c('sf','sfc'))) {
    poly <- st_geometry(poly) # extract sfc if not already sfc
    geomtype <- st_geometry_type(poly, by_geometry = FALSE)
    if (geomtype == 'GEOMETRY') {   # 2023-06-02
        geomtype <- st_geometry_type(poly, by_geometry = TRUE)
    }
    if (!all(geomtype %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop ("poly sf/sfc should be of type POLYGON or MULTIPOLYGON")
    }
    poly
  }
  else if (inherits(poly, 'SpatialPolygons')) {   # also SPDF?
    st_as_sfc(poly)
  }
  else if (inherits(poly, 'SpatVector')) {
    st_as_sfc(as(poly,"Spatial"))
  }
  else if (inherits(poly, c('matrix', 'data.frame'))) {
    ## input is 2-column matrix for a single polygon
    poly <- matrix(unlist(poly), ncol = 2)
    poly <- rbind (poly, poly[1,])  ## force closure of polygon
    st_sfc(st_polygon(list(poly)))
  }
  else stop (class(poly), " not valid input to boundarytoSF")
}

#-------------------------------------------------------------------------------

pointsInPolygon <- function (xy, poly, logical = TRUE) {
  # xy is 2-column matrix or data.frame of coordinates
  if (inherits(poly, 'mask')) { 
    if (ms(poly))
      stop ("multi-session masks not supported")
    sp <- spacing(poly)
    minx <- min(poly$x, na.rm = TRUE)
    miny <- min(poly$y, na.rm = TRUE)
    mask <- sweep(poly, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
    mask <- round(mask/sp) + 1
    xy <- matrix(unlist(xy), ncol = 2)  ## in case dataframe
    xy <- sweep(xy, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
    xy <- round(xy/sp) + 1
    xy[xy<=0] <- NA
    xy[,1][xy[,1]>max(mask$x, na.rm = TRUE)] <- NA
    xy[,2][xy[,2]>max(mask$y, na.rm = TRUE)] <- NA
    
    maskmatrix <- matrix(0, ncol = max(mask$y, na.rm = TRUE), nrow = max(mask$x, na.rm = TRUE))
    maskmatrix[as.matrix(mask)] <- 1:nrow(mask)
    inside <- maskmatrix[as.matrix(xy)]
    inside[is.na(inside)] <- 0
    if (logical)
      inside <- inside > 0
    inside
  }
  else {
    poly <- boundarytoSF(poly)
    if (inherits(poly, c('sf','sfc'))) {
      xy <- st_as_sf(data.frame(xy), coords = 1:2)
      st_crs(xy) <- st_crs(poly)
      apply(st_within(xy, poly, sparse = FALSE), 1, any)
    }
    else {
      stop ("unknown input to pointsInPolygon")
    }
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

maskboolean <- function (ch, mask, threshold) {
  if (ms(ch)) {
    if (!ms(mask)) stop ("masklookup: multisession ch requires multisession mask")
    outlist <- mapply(maskboolean, ch, mask, MoreArgs = list(threshold = threshold), SIMPLIFY = FALSE)
    outlist
  }
  else {
    id <- animalID(ch, names = FALSE, sortorder = 'snk')
    tr <- trap(ch, names = FALSE, sortorder = 'snk')
    trps <- traps(ch)
    m <- nrow(mask)
    if (!is.null(threshold) && all(detector(trps) %in% .localstuff$pointdetectors)) {
      df <- data.frame(id = id, x = trps$x[tr], y = trps$y[tr])
      x <- tapply(df$x, df$id, mean, na.rm=T)
      y <- tapply(df$y, df$id, mean, na.rm=T)
      xy <- data.frame(x=x,y=y)
      d2 <- edist(xy, mask)
      out <- (d2 <= threshold^2)
    }
    else {
      ## NULL option
      out <- matrix(TRUE, nrow = nrow(ch), ncol = m)
    }
    out
  }
}

#-------------------------------------------------------------------------------

uniquerownames <- function (capthist) {
    if (!ms(capthist)) {
        return(capthist)
    }
    else {
        last <- 0
        for (i in 1:length(capthist)) {
            nr <- nrow(capthist[[i]])
            if (nr > 0) {
            rownames(capthist[[i]]) <- last + (1:nr)
            last <- last+nr
            }
        }
        capthist
    }
}

#-------------------------------------------------------------------------------

selectCHsession <- function(capthist, sessnum) {
    if (ms(capthist)) 
        capthist[[sessnum]]
    else 
        capthist
}

#-------------------------------------------------------------------------------

stringsAsFactors <- function (DF) {
    # convert any character columns of a data.frame (or list) to factor
    if (is.list(DF) && length(DF)>0) {    ## bug fix 2020-08-14
        chr <- sapply(DF, is.character)
        DF[chr] <- lapply(DF[chr], as.factor)
    }
    DF
}

#-------------------------------------------------------------------------------

## function to assign all-ones usage matrix
uniformusage <- function(object, noccasions) {
  if (inherits(object, 'capthist')) {
    if (ms(object)) {
      for (r in 1:length(object)) {
        ndet <- dim(object[[r]])[3]
        noccasions <- dim(object[[r]])[2]
        usage(traps(object[[r]])) <- matrix(1, ndet, noccasions)
      }
    }
    else {
      ndet <- dim(object)[3]
      noccasions <- dim(object)[2]
      usage(traps(object)) <- matrix(1, ndet, noccasions)
    }
  }
  else if (inherits(object, 'traps')) {
    if (missing(noccasions)) {
      stop ('noccasions should be specified for traps input')
    }
    if (ms(object)) {
      for (r in 1:length(object)) {
        ndet <- ndetector(object[[r]])
        usage(object[[r]]) <- matrix(1, ndet, noccasions)
      }
    }
    else {
      ndet <- ndetector(object)
      usage(object) <- matrix(1, ndet, noccasions)
    }
  }
  object
}

#-------------------------------------------------------------------------------

sfrotate <- function (x, degrees, centrexy = NULL, usecentroid = FALSE) {
    rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
    gx <- st_geometry(x)
    if (is.null(centrexy)) {
        if (usecentroid) {
            centrexy <- st_centroid(gx)[1,]   # unique centre
        }
        else {
            centrexy <- st_centroid(st_as_sfc(st_bbox(x)))
        }
    } 
    else {
        centrexy <- st_sfc(st_point(centrexy) )
    }
    (gx - centrexy) * rot(degrees/360*2*pi) + centrexy
}

#-------------------------------------------------------------------------------

# Based on Tim Salabim stackoverflow Jul 12 2018
# https://stackoverflow.com/questions/51292952/snap-a-point-to-the-closest-point-on-a-line-segment-using-sf

snap_points <- function(x, y, max_dist = 1000) {
    
    if (inherits(x, "sf")) n = nrow(x)
    if (inherits(x, "sfc")) n = length(x)
    
    out = do.call(c,
        lapply(seq(n), function(i) {
            nrst = st_nearest_points(st_geometry(x)[i], y)
            nrst_len = st_length(nrst)
            nrst_mn = which.min(nrst_len)
            if (as.vector(nrst_len[nrst_mn]) > max_dist) return(st_geometry(x)[i])
            return(st_cast(nrst[nrst_mn], "POINT")[2])
        })
    )
    return(out)
}

#-------------------------------------------------------------------------------

im2mask <- function(im) {
    # spatstat im object to mask
    df <- as.data.frame(im)
    names(df) <- c('x','y','Lambda')
    df$Lambda <- df$Lambda * 1e4   # per hectare
    read.mask(data = df, spacing = im$xstep)
}
#-------------------------------------------------------------------------------
