library(rstrauss)
set.seed(1)
d3 <- F
comp <- 0.5
bb <- cbind(c(0,comp) -comp/2, c(0,1/comp)-0.5/comp)
bbm <- cbind(0:1, 0:1)-0.5
M <- diag(c(1/comp, comp))

if(d3) bb<-cbind(bb, 0:1)
pp0 <- rstrauss(500, .01, R<-0.03, perfect=TRUE, bbox=bb, iter=2e4)
xm <- pp0$x%*%M
pp <- list(x=xm, bbox=bbm)
