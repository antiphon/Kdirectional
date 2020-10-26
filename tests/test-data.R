library(rstrauss)
set.seed(1)
if(!exists("d3"))d3 <- F
comp <- .8
M <- if(d3) diag(c(1/sqrt(comp), comp, comp)) else diag(c(1/comp, comp))




bb <- cbind(0:1,0:1)*2-1
if(d3) bb <- cbind(bb, 0:1)
pp0 <- rstrauss(500, .01, R<-0.09, perfect=TRUE, bbox=bb, iter=2e4)
xm <- coord_affine(pp0$x, M)
bbm <- bbox_affine(bb, M)
pp <- list(x=xm, bbox=bbm)

