#' test histogram limits
library(rstrauss)
library(Kdirectional)
beta <- 10000
R <- strauss_theory(beta, 2)*0.8
#set.seed(1)
x0 <- rstrauss(beta, .015, R, bbox=cbind(0:1,0:1), perfect=T, iter = 1e5)

library(icetools)
p1 <- function(x, ... , main=""){
  x$domain <- list(bbox=x$bbox)
  x <- deform.pp3d(x, ...)
  angles <- nnangle(x$x)[which(bbox.distance(x$bbox, x$x)>0.1),]
  r<-angles.flower(angles[,1], ci=T, k = k<-25)
  cdf <- angles.cdf(angles[,1], ci=T)
  plot(cdf$ang, cdf$cdf-(v<-cdf$ang/(2*pi)), "l", ylim=yl<-c(-1,1)*0.1, main=main)
  with(cdf, {lines(angle, CI5-v, col=2);lines(angle,CI95-v, col=2)})
  abline(v=seq(0,2*pi, by=pi/4), col="gray90")
}
#' compress

par(mfrow=c(4, 2))
p1(x0, dy=1, main="no compression")
p1(x0, dy=0.9, main="dy=0.9")
p1(x0, dy=0.7, main="dy=0.7")
p1(x0, dy=0.5, main="dy=0.5")
