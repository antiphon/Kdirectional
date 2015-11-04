# check the v2

library(devtools)
load_all(".")


if(!exists("x1")){
  library(rstrauss)
  w <- c(-1,1)*.5
  bb <- cbind(w,w)
  comp <- 0.7
  bbsim <- cbind(w*comp,w/comp)
  z <- rstrauss(250, 0.1, 0.04, bbox = bbsim, verb=T, perfect=T, iter=1e5)
  y <- z$x %*%diag(c(comp,1/comp))
  x1 <- list(x=y, bbox=bb)
}

n <- 10

el <- fry_ellipsoids(x1, keep_data = T, nvec =1:5, double=F, nangles=n)
ell <- fry_ellipsoids(x1, keep_data = T, nvec =1:5, double=T, nangles=n)
el2 <- fry_ellipsoids2(x1, keep_data = T, nvec =1:5, double=F, nangles=n)
ell2 <- fry_ellipsoids2(x1, keep_data = T, nvec =1:5, double=T, nangles=n)

par(mfrow=c(2,2))
pp <- function(e) plot(e, zoom=.2)#, sectors=c(1:nrow(el$grid_unit)))
pp(el)
pp(ell)
pp(el2)
pp(ell2)
