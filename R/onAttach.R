###############################################################################
## package 'secrpoly'
## onAttach.R
## last changed 2024-01-29
###############################################################################

.onAttach <- function (libname, pkgname) {
    version <- paste0(packageVersion('secrpoly'), .localstuff$packageType)
    packageStartupMessage( "This is secrpoly ", version,
                           ". For overview type ?secrpoly" )
}

## .onLoad is preferred if actions are required for single functions 
## that may be called without attaching package