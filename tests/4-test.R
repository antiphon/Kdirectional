#' test
library(spatstat)
library(Kdirectional)

p <- rStrauss(1000, 0.1, 0.03, W=owin(c(0,1),c(0,2)))
k <- Fest(p, r=seq(0, 0.1, length=10))
plot(k, cbind(theo,rs)~r, xlim=c(0,.1))
#'
pp <- list(x=cbind(p$x, p$y), bbox=cbind(c(0,1),c(0,2)))
k2 <- Fest_directional(pp, u=c(1,0), theta=pi, r=k$r)
lines(k2$r, k2$rs, col=3)
#'

