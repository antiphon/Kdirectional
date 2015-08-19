#' anisotropic pcf
library(devtools)
load_all(".")
library(rgl)
library(icetools)
library(rstrauss)
library(plotrix)
set.seed(2)
d3 <- F
bb <- bbox.default(0:1, 0:1, if(d3) 0:1 else NULL)-0.5
co <- 0.7
bbm <- deform.bbox(bb, dz=1/co)

pp <- rstrauss(1000, .01, R<-0.01, perfect=TRUE, bbox=bbm)
#ppm <- list(x=deform.coords(pp$x, dy=co), bbox=bb)
#'
dirs <- check_directions(n_dir=15, dim = 2)

g <- pcf_anisotropic(pp, theta=dirs$theta, r=seq(0, 0.1, length=50))
#gm <- pcf_anisotropic(ppm, h=g$h, r=g$r, theta=g$theta)

if(ncol(bb)==2){
  g$epsilon <- dirs$epsilon
  class(g) <- c("pcf_sector", is(g))
  sectorplot(g, overlapfactor = 0.38, zlim=c(0,2))
#plot(gm)
}
