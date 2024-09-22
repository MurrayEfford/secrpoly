## 2024-09-23

library(secrpoly)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

####################################################################

# Polygon detectors

set.seed(123)
polyX <- make.poly (exclusive = TRUE)
poly  <- make.poly (exclusive = FALSE)

CHpolyX <- sim.capthist (polyX, popn = list(D = 100, buffer = 100), 
                        detectfn = 'HHN', detectpar = list(lambda0 = 0.5, sigma = 25), 
                        noccasions = 1)

CHpoly <- sim.capthist (poly, popn = list(D = 100, buffer = 100), 
                        detectfn = 'HHN', detectpar = list(lambda0 = 0.5, sigma = 25), 
                        noccasions = 1)

test_that("correct likelihood (polygonX data)", {
    args <- list(capthist = CHpolyX, detectfn = 'HHN', buffer = 100, 
                 start = log(c(100, 0.5,  25)), verify = FALSE,
                 details = list(LLonly = TRUE))
    LL1 <- do.call(secrpoly.fit, args)[1]
    expect_equal(LL1, -510.65755, tolerance = 1e-4, check.attributes = FALSE)
    
})

test_that("correct likelihood (polygon data)", {
    args <- list(capthist = CHpoly, detectfn = 'HHN', buffer = 100, 
                 start = log(c(100, 0.5,  25)), verify = FALSE,
                 details = list(LLonly = TRUE))
    LL1 <- do.call(secrpoly.fit, args)[1]
    expect_equal(LL1, -482.73604, tolerance = 1e-4, check.attributes = FALSE)
    
})

###############################################################################
