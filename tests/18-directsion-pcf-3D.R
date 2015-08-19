#' anisotropic pcf 3d

library(devtools)
load_all(".")
library(rstrauss)
library(rgl)

set.seed(1)
comp <- .7
comps <- sqrt(comp)
bb <- cbind(c(0,comps) -comps/2, c(0,comps) -comps/2, c(0,1/comp)-0.5/comp)
bbm <- cbind(0:1, 0:1, 0:1)-0.5
M <- diag(c(1/comps, 1/comps, comp))

pp0 <- rstrauss(500, .01, R<-0.1, perfect=TRUE, bbox=bb, iter=2e4)

xm <- pp0$x%*%M

pp <- list(x=xm, bbox=bbm)


g <- pcf_directions(pp0, n_dir=2, r=seq(0,0.2, length=10))

co <- values2colors(g$est, col=rainbow, zlim=c(0,2))

plot3d(g$directions, col=co, aspect=F, size=.1)
spheres3d(g$directions, col=co, radius=rr<-0.005)
spheres3d(-g$directions, col=co, radius=rr)

play3d(spin3d())
