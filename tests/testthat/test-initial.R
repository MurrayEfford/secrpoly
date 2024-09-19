## Started 2024-09-19
library(secrpoly)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################
## Multi-polygon bug secr 2021-05-18

datadir <- system.file("extdata", package = "secrpoly")
polyexample1 <- read.traps(file = paste0(datadir, '/polygonexample1.txt'), 
    detector = 'polygon')
polygonCH <- sim.capthist(polyexample1, popn = list(D = 1, buffer = 200),
    detectfn = 'HHN', detectpar = list(lambda0 = 5, sigma = 50),
    noccasions = 1, seed = 123)

test_that("simulated polygon data have correct RPSV", {
    rpsv <- RPSV(polygonCH, CC = TRUE)
    expect_equal(rpsv, 45.16401, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct likelihood (multi-detector polygon data)", {
    args <- list(capthist = polygonCH, buffer = 200, detectfn = 'HHN',
        start = list(D=1, lambda0=5, sigma = 50), verify = FALSE, 
        details = list(LLonly = TRUE))
    LL <- do.call(secrpoly.fit, args)[1]
    expect_equal(LL, -2026.92169, tolerance = 1e-4, check.attributes = FALSE)
})
###############################################################################
