#############################################################################
## package 'secrpoly'
## Methods for classes traps, capthist and mask
## 2024-01-29
###############################################################################

## traps object
polyID <- function (object)    {
    if (ms(object)) {
        polyID (object[[1]])
    }
    else {
        if (inherits(object,'traps')) {
            temp <- attr(object,'polyID')
            if (is.null(temp)) temp <- factor(1:nrow(object))   ## all different
            temp
        }
        else
        if (inherits(object,'capthist')) {
            stop ("use trap() to extract polyID from 'capthist' object")
        }
        else stop ("polyID requires 'traps' object")
    }
}

## traps object
transectID <- function (object)    {
    if (ms(object)) {
        transectID (object[[1]])
    }
    else {
        if (inherits(object,'traps')) {
            if (!all(detector(object) %in% c('transect','transectX')))
                stop ("requires transect detector")
            temp <- attr(object,'polyID',exact = TRUE)
            if (is.null(temp)) temp <- factor(1:nrow(object))
            temp
        }
        else
        if (inherits(object,'capthist')) {
            stop ("use trap() to extract transectID from 'capthist' object")
        }
        else stop ("transectID requires 'traps' object")
    }
}

xy <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, xy)
    }
    else {
        if (any(detector(traps(object)) %in%
            c('polygonX', 'transectX', 'polygon','transect'))) {
            attr(object, 'detectedXY',exact = TRUE)
        }
        else
            NULL
    }
}

alongtransect <- function (object, tol = 0.01) {
    ptalongtransect <- function (i) {
        ## where is point i on its transect k?
        k <- trans[i]
        transectxy <- as.matrix(lxy[[k]])
        nr <- nrow(transectxy)
        alongtransectcpp (
            as.matrix (xyi[i,,drop=FALSE]),
            as.matrix (transectxy),
            as.integer (0),
            as.integer (nr-1),
            as.double (tol))
    }
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    
    if (ms(object)) {
        lapply(object, alongtransect, tol = tol)
    }
    else {
        trps <- traps(object)
        
        if (all(detector(trps) %in% c('transectX', 'transect'))) {
            
            #             trans <- trap(object, names = TRUE)
            #             xyi <- xy(object)
            #             ## 2015-09-02 change to fix occsim: remove 'S' prefix
            #             lxy <- split (trps, levels(transectID(trps)), prefix = "")
            
            trans <- trap(object, names = FALSE, sortorder = 'ksn')
            xyi <- xy(object)
            lxy <- split (trps, transectID(trps))
            
            sapply(1:nrow(xyi), ptalongtransect)
        }
        else
            NULL
    }
}

# under development 2021-12-11, 2022-02-01
distancetotransect <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    
    if (ms(object)) {
        lapply(object, distancetotransect)
    }
    else {
        trps <- traps(object)
        
        if (all(detector(trps) %in% c('transectX', 'transect'))) {
            
            xyi <- xy(object)
            
            # split by transect
            lxy <- split (trps, transectID(trps))
            vlist <- lapply(lxy, as.matrix)
            
            # each transect as sfg
            vlist <- lapply(vlist, st_linestring)
            # combine linestrings in one sfc
            v <- st_sfc(vlist)
            xy <- st_as_sf(as.data.frame(xyi), coords = 1:2)
            neari <- st_nearest_feature(xy, v)
            st_distance(xy,v)[,neari]   # distance to nearest feature

        }
        else
            NULL
    }
}

searcharea <- function (object)    {
    ## requires traps object
    ## discarded some obsolete code 2016-10-08
    if (ms(object)) {
        ## 2011-06-24
        lapply(object, searcharea)
    }
    else {
        if (all(detector(object) %in% c('polygon','polygonX'))) {
            ## assume poly closed
            sapply(split(object, levels(polyID(object))), polyarea)
        }
        else NA
    }
}

transectlength <- function (object)    {
## requires traps object
    if (ms(object)) {
        ## 2011-06-24
        lapply(object, transectlength)
    }
    else {
        if (all(detector(object) %in% c('transect','transectX'))) {
            calclength <- function (xy) {
                nr <- nrow(xy)    ## number of vertices
                segments <- (xy$x[-nr] - xy$x[-1])^2 + (xy$y[-nr] - xy$y[-1])^2
                sum (segments^0.5)
            }
            sapply(split(object, polyID(object)), calclength)
        }
        else NA
    }
}


'polyID<-' <- function (object, value) {
    if (!inherits(object,'traps'))
        warning ("polyID requires 'traps' object")
    if (length(value)==1) value <- rep(value, nrow(object))
    value <- factor(value)
    structure (object, polyID = value)
}

'transectID<-' <- function (object, value) {
    if (length(value)==1) value <- rep(value, nrow(object))
    value <- factor(value)
    structure (object, polyID = value)
}

'xy<-' <- function (object, value) {
    if (!is.null(value)) {
        polyoccasions <- expanddet(object) %in% .localstuff$polydetectors
        ndetections <- sum(abs(object[,polyoccasions,]))
        if (nrow(value) != ndetections)
            stop ("requires one location per detection")
        if ((sum(polyoccasions)==0) |
            !(inherits(object,'capthist')))
            stop ("requires 'capthist' object with ",
                  "polygon-like detector")
        if (ms(object))
            stop ("requires single-session 'capthist' object")
    }
    structure (object, detectedXY = value)
}
