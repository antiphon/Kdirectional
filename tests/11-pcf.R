#' test pcf dir

library(devtools)
load_all(".")
library(rstrauss)
library(spatstat)
set.seed(1)

p <- rStrauss(100, 0.1, 0.03)
k <- pcf(p, r=seq(0, 0.2, length=50))
#'
pp <- list(x=cbind(p$x, p$y), bbox=cbind(c(0,1),c(0,1)))

k2 <- pcf_sector(pp,  r=k$r, theta=(pi), epsilon=pi)

k3 <- pcf_sector(pp, r=k$r, theta=k2$theta, epsilon=k2$epsilon, correction="none")


plot(k, trans~r, xlim=c(0,.1), main="pcf", lwd=4)
lines(k2$r, k2$theo, col=4)
lines(k2$r, k2$est[,1], col=2, lty=2)
lines(k2$r, k3$est[,1], col=3, lwd=3, lty=2)

sectorplot(k2, zlim=c(0,2), res = 12)
