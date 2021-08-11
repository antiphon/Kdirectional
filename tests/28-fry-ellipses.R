# fry ellipsoid finding
library(devtools)
load_all(".")

library(rstrauss)
library(ellipsoid)

if(!exists("x1")){
  w <- c(-1,1)*.5
  bb <- cbind(w,w)
  comp <- 0.5
  bbsim <- cbind(w*comp,w/comp)
  z <- rstrauss(100, 0.001, 0.08, bbox = bbsim)
  y <- z$x %*% solve(diag(c(comp,1/comp)))
  x1 <- list(x=y, bbox=bb)
  x1 <- spatstat.geom::ppp(y[,1], y[,2], spatstat.geom::as.owin(c(bb)))
}

el <- fry_ellipsoids(x1, nvec =2:7)
ave <- mean_ellipsoids(el$ellipsoid, nsim = 10000)

# cylindric instead of sector:
el2 <- fry_ellipsoids(x1, nvec =2:7, cylindric=T, eps=.04)
ave2 <- mean_ellipsoids(el2$ellipsoids, nsim=2000)


par(mfrow=c(2,1))
plot.fryellipsoids(el, T, T, zoom = .2)
plot.ellipsoid(ave, lwd=5)
plot.fryellipsoids(el2, T, zoom = .5, sectors = c(2, 7))
plot.ellipsoid(ave2, lwd=5)
is(ave)
