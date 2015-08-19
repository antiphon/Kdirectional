#' test
library(spatstat)
library(Kdirectional)

p <- rStrauss(100, 0.1, 0.03)
k <- Kest(p, r=seq(0, 0.2, length=50))
plot(k, trans~r, xlim=c(0,.1))

#'
pp <- list(x=cbind(p$x, p$y), bbox=cbind(c(0,1),c(0,1)))
k2 <- Kest_directional(pp, u=c(0, 1), theta=pi, r=k$r)
lines(k2$r, k2$trans, col=2)
lines(k2$r, k2$theo, col=3)
