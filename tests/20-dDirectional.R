#' Test directional distribution
#' 

library(devtools)
load_all(".")
library(sphere)
library(rstrauss)
bb <- cbind(0:1, 0:1)
x <- rstrauss(10000, 0.1, 0.005, bbox=bb)

z <- 0.5
M <- diag(c(1/z, z))

x2 <- list(x=x$x%*%M)
x2$bbox <- apply(x2$x, 2, range)
#' sanity
g0 <- dDirectional(x, n_dir=10, r=seq(0,0.08, length=100))
g1 <- dDirectional(x2, n_dir=10, r=seq(0,0.08, length=100))


par(mfrow=c(1,1))
plot(g0$r, g0$est[,1], "l", ylim=c(0,1))
i<-0
for(gd in list(g0, g1)){
apply(gd$est, 2, lines, x=gd$r, col=i<-i+1)
}

theo <-  1-exp(-nrow(x$x)*g0$epsilon*g0$r^2)
lines(g0$r, theo, lwd=3, col=3)
