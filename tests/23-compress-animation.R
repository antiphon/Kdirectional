#' test compress function

library(devtools)
load_all(".")


library(rstrauss)
z <- 0.6
bb0 <- cbind(c(0,z), c(0,1/z))

bb1 <- cbind(c(0,z), c(0,1/z))

p <- rstrauss(100, 0.01, 0.1, bbox=bb0)

x <- p$x
for(i in seq(0,2, length=15)){
  ang <- pi/2
  
  u <- (1 - 0.4*(i/2))*c(cos(ang), sin(ang))
  
  xx <- compress(x, u, vol_preserve = F, F)
  
  
  par(mfrow=c(2,2))
  cr <- function(){abline(h=0, col="gray60");abline(v=0, col="gray60")}
  plot(x, asp=1, col=3, xlim=Lx<-c(0,1), ylim=Ly<-c(0,1), pch=19, cex=ce<-0.4)
  cr()
  points(u[1],u[2], col=2)
  
  plot(xx, asp=1, col=3, xlim=Lx, ylim=Ly, pch=19, cex=ce)
  points(u[1],u[2], col=2)
  cr()
  plot(x, asp=1, col=3, xlim=Lx, ylim=Ly, pch=19, cex=ce)
  cr()
  points(xx, col=4, pch=19, cex=.2)
  for(i in 1:nrow(x)){
    arrows(x[i,1], x[i,2], xx[i,1] ,xx[i,2], col="gray50", length = 0.1)
  }
}