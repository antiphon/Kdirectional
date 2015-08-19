# fry ellipsoid finding
library(devtools)
load_all(".")

library(rstrauss)
library(sphere)
if(!exists("x1")){
  w <- c(-1,1)*.5
  bb <- cbind(w,w)
  comp <- 0.5
  bbsim <- cbind(w*comp,w/comp)
  z <- rstrauss(100, 0.001, 0.08, bbox = bbsim)
  y <- compress(z$x, u=c(0,1/comp))
  x1 <- list(x=y, bbox=bb)
}

el <- fry_ellipsoids(x1, nvec =2:7)
ave <- mean_ellipsoids(el$ellipsoids)

# cylindric instead of sector:
el2 <- fry_ellipsoids(x1, nvec =2:7, cylindric=T, eps=.04)
ave2 <- mean_ellipsoids(el2$ellipsoids, nsim=2000)


par(mfrow=c(2,1))
plot(el, T, T, sector=F)
plot(ave, lwd=5)
plot(el2, T, sector=c(2,9))
plot(ave2, lwd=5)