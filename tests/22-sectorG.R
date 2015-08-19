#' Test sector-G
#' 

#' Test directional distribution
#' 

library(devtools)
load_all(".")
library(sphere)
library(rstrauss)
library(plotrix)
bb <- cbind(0:1, 0:1)
x <- rstrauss(1000, 0.1, 0.01, bbox=bb)

z <- 0.6
M <- diag(c(1/z, z))

x2 <- list(x=x$x%*%M)
x2$bbox <- apply(x2$x, 2, range)
#' sanity
g0 <- Gsector(x, n_dir=12, r=seq(0,.2, length=200))
g1 <- Gsector(x2, theta=g0$theta, epsilon=g0$epsilon, r=g0$r)


par(mfrow=c(2,2))
 

i<-0
for(g in list(g0, g1)){
  i<-i+1
  theo <- g$theo
  plot(g$r, g$est[,1]-theo, "l", ylim=c(-1,1)*0.05)
  apply(g$est, 2, function(v) lines(v-theo, x=g$r, col=i))
  abline(h=0, col=3)
}
abline(h=0, col=3)

sectorplot(g0)
sectorplot(g1)



