#' anisotropic inhomogeneous K, cylinder version

library(devtools)
load_all(".")
library(rstrauss)

if(!exists("pp")){
  set.seed(1)
  dim <- 2
  comp <- 0.7
  side <- c(0,1)
  bbm <- matrix(rep(side, dim), nc=dim)
  M <- diag(c(rep(1/comp^(1/(dim-1)), dim-1), comp))
  bb <- bbm%*%solve(M)
  pp0 <- rstrauss(n=dim*100, gamma = .01, range = R<-0.04, bbox=bb, iter=1e5, verb=T, toroidal=T)
  xm <- pp0$x%*%M
  pp <- list(x=xm, bbox=bbm)
}

# sanity checks
if(0){
  l <- rep(nrow(pp$x)/bbox_volume(pp$bbox), nrow(pp$x))
  
  p <- Kest_anin_cylinder(pp, lambda=l, epsilon = e<-0.03)
  pd <- Kest_anin_cylinder(pp, lambda=l, epsilon = e*2)
  
  # compare to sector
  p2 <- pcf_anin(pp, lambda=l, epsilon = e2<-pi/8)
  pd2 <- pcf_anin(pp, lambda=l, epsilon = e2*2)
  
  par(mfrow=c(2+1*(dim==2),2))
  plot(p, main=e)
  plot(pd, main =e*2)
  plot(p2, main=e2)
  plot(pd2,  main =e2*2)
  if(dim==2) plot(pp$x, asp=1)
}

# adaptive
if(1) {
  l <- rep(nrow(pp$x)/bbox_volume(pp$bbox), nrow(pp$x))
  
  p  <- Kest_anin(pp, lambda=l)
  pa <- Kest_anin_cylinder(pp, lambda=l, aspect = 1/3)
  
  par(mfrow=c(2,1))
  plot(p, log = "y")
  plot(pa, log = "y")
}