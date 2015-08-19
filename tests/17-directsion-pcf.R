#' anisotropic pcf

library(devtools)
load_all(".")
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
#plot(pp0$x, asp=1); points(pp$x, col=2)
#'

g0 <- pcf_directions(pp0, n_dir=25, r=seq(0,0.1, length=30))
g <- pcf_directions(pp, n_dir=25, r=seq(0,0.1, length=30))

#' compare:
par(mfrow=c(1,2))
flower.pcf_directions(g0)
flower(g)

