## code used in debugging to source functions into workspace
## requires previous modification of NAMESPACE to exportPattern("*cpp")
## 2023-10-17 moved from secr4 folder

library(secrpoly)
library(RcppParallel)
library(RcppNumerical)

source('d:/density secr 4.6/secrpoly/R/RcppExports.R')
source('d:/density secr 4.6/secrpoly/R/secrpoly.fit.R')
source('d:/density secr 4.6/secrpoly/R/utility.R')
source('d:/density secr 4.6/secrpoly/R/generalsecrloglik.R')
source('d:/density secr 4.6/secrpoly/R/loglikhelperfn.R')
source('d:/density secr 4.6/secrpoly/R/reparameterize.R')
source('d:/density secr 4.6/secrpoly/R/preparedata.R')
source('d:/density secr 4.6/secrpoly/R/fxi.R')
