# testing.R

library(secrpoly)
library(secr)
# source('d:/density secr 4.6/secrpoly/sourcecode.R')

fit0 <- secr.fit(hornedlizardCH, buffer=100, ncores=18)

plot(hornedlizardCH, tracks=T, border=30)
fxi.contour(fit0, 1:68, add=T)
