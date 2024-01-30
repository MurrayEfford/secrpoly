## Started 2024-01-29

library(secrpoly)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

####################################################################

# Poisson counts

set.seed(123)
detectors <- make.grid (nx = 6, ny = 8, detector = "count")
CHpois <- sim.capthist (detectors, popn = list(D = 10, buffer = 100), 
    detectpar = list(g0 = 0.2, sigma = 25), noccasions = 1)

test_that("correct likelihood (Poisson count data)", {
    args <- list(capthist = CHpois, detectfn = 'HN', buffer = 100, 
        start = c(2.121835788, -1.040097594,  3.201521728), verify = FALSE)
    
    args$details <- list(LLonly = TRUE, fastproximity = TRUE)
    LL1 <- do.call(secr.fit, args)[1]
    expect_equal(LL1, -122.138538, tolerance = 1e-4, check.attributes = FALSE)
    
})

###############################################################################
