#' anisotropic inhomogeneous pcf, cylinder version

library(devtools)
load_all(".")
library(rstrauss)

if(!exists("pp")){
  set.seed(1)
  dim <- 3
  comp <- 0.7
  side <- c(1,2)
  bbm <- matrix(rep(side, dim), nc=dim)
  M <- diag(c(rep(1/comp^(1/(dim-1)), dim-1), comp))
  bb <- bbm%*%solve(M)
  pp0 <- rstrauss(n=dim*100, gamma = .01, range = R<-0.03, bbox=bb, iter=1e5, verb=T, toroidal=T)
  xm <- pp0$x%*%M
  pp <- list(x=xm, bbox=bbm)
}
l <- rep(nrow(pp$x)/bbox_volume(pp$bbox), nrow(pp$x))
p <- pcf_anin_cylinder(pp, lambda=l, epsilon = e<-0.05, stoyan=f<-0.3)
pd <- pcf_anin_cylinder(pp, lambda=l, epsilon = e*2, stoyan=f)

# compare to sector
p2 <- pcf_anin(pp, lambda=l, epsilon = e2<-pi/8,  divisor="d",stoyan=f<-0.3)
pd2 <- pcf_anin(pp, lambda=l, epsilon = e2*2,  divisor="d",stoyan=f)

par(mfrow=c(2+1*(dim==2),2))
plot(p, main=e)
plot(pd, main =e*2)
plot(p2, main=e2)
plot(pd2, main =e2*2)
if(dim==2) plot(pp$x, asp=1)
