# test the wrapper function pcf_directional

library(devtools)
load_all(".")
#library(rstrauss)
library(spatstat)
set.seed(1)

p <- rStrauss(100, 0.1, 0.03)
k <- pcf(p, r=seq(0, 0.2, length=50))
#'
pp <- list(x=cbind(p$x, p$y), bbox=cbind(c(0,1),c(0,1)))

k2 <- pcf_directional(pp)
