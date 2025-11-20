################################################################################
## package 'secrpoly'
## utility.R
## 2024-01-29, 2024-09-20
## 2025-11-21 removed duplicates of secr utility functions
################################################################################

# Global variables in namespace
#
## define a local environment for temporary variables e.g. iter
.localstuff <- new.env()

.localstuff$packageType <- ' pre-release'
##.localstuff$packageType <- ''

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

.localstuff$learnedresponses <- c('b', 'bk', 'B', 'k', 'Bk') 

#-------------------------------------------------------------------------------

secrpoly_valid.pnames <- function (details, CL, detectfn, alltelem, sighting, nmix) {
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

memo <- function (text, trace) {
    ## could use message(text), but does not immediately flush console
    if (trace) { cat (text, '\n')
    flush.console() }
}
#-------------------------------------------------------------------------------

var.in.model <- function(v,m) v %in% unlist(lapply(m, all.vars))

#-------------------------------------------------------------------------------
