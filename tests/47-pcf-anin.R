#' anisotropic inhomogeneous pcf, sector version

library(devtools)
load_all(".")
library(rstrauss)

if(!exists("pp")){
  set.seed(1)
  d3 <- F
  comp <- 0.9
  #bb <- cbind(c(0,comp) -comp/2, c(0,1/comp)-0.5/comp)
  side <- c(1,2)
  side2<- c(1,3)
  bbm <- cbind(side2,side)
  M <- diag(c(1/comp, comp))
  if(d3)M <- diag(c(sqrt(1/comp),sqrt(1/comp), comp))
  if(d3) bbm<-cbind(bbm, side)
  bb <- bbm%*%solve(M)
  pp0 <- rstrauss(500, .01, R<-0.03+ifelse(d3,0.08,0), perfect=TRUE, bbox=bb, iter=2e4)
  xm <- pp0$x%*%M
  pp <- list(x=xm, bbox=bbm)
}
l <- rep(nrow(pp$x)/bbox_volume(pp$bbox), nrow(pp$x))
p <- pcf_anin(pp, lambda=l, epsilon = e<-pi/4, divisor = "d", stoyan=f<-0.3)
pd <- pcf_anin(pp, lambda=l, epsilon = e, divisor="r", stoyan=f)

par(mfrow=c(1,2))
plot(p, main="divisor d")
plot(pd, main =" divisor r")

