#' anisotropic pcf as by old school boys

library(devtools)
load_all(".")
library(rstrauss)

if(!exists("pp")){
  set.seed(1)
  d3 <- T
  comp <- 0.6
  #bb <- cbind(c(0,comp) -comp/2, c(0,1/comp)-0.5/comp)
  side <- c(0,1)
  side2<- c(0,3)
  bbm <- cbind(side2,side)
  M <- diag(c(1/comp, comp))
  if(d3)M <- diag(c(sqrt(1/comp),sqrt(1/comp), comp))
  if(d3) bbm<-cbind(bbm, side)
  bb <- bbm%*%solve(M)
  pp0 <- rstrauss(500, .01, R<-0.03+ifelse(d3,0.08,0), perfect=TRUE, bbox=bb, iter=2e4)
  xm <- pp0$x%*%M
  pp <- list(x=xm, bbox=bbm)
}


r = seq(0,.9, l = 100)
l <- rep(nrow(pp$x)/bbox_volume(pp$bbox), nrow(pp$x))
d <- ncol(pp$x)
#
pd <- pcf_anisotropic(pp, r = r, h = c(0.35/l[1]^(1/d), e <- pi/5))

# compare
p <- pcf_anin(pp, lambda=l, r=r, epsilon = e, stoyan=f<-0.35)
pc <- pcf_anin_cylinder(pp, lambda=l, r=r, epsilon = 0.08, stoyan=f)

par(mfrow=c(2,2))
plot(p, main="sector")
plot(pc, main="cylinder")
plot(pd, main ="anisotropic polar")

