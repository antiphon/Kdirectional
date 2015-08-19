#' test compress function

library(devtools)
load_all(".")

d3 <- TRUE
# 2D #####
if(!d3){
library(rstrauss)
z <- 0.9
bb0 <- cbind(c(-1,1), c(-1,1))

p <- rstrauss(100, 0.01, 0.1, bbox=bb0)

cr <- function(){abline(h=0, col="gray60");abline(v=0, col="gray60")} # origin lines

x <- p$x
par(mfrow=c(5,3))
angs <- seq(0, pi, length=15)
for(a in angs){
  ang <- a
  #' force direction
  u <- z*c(cos(ang), sin(ang))
  
  xx <- compress(x, u, vol_preserve = T, F)
  
  plot(xx, asp=1, col=3, xlim=Lx<-c(-2,2), ylim=Ly<-c(-2,2), pch=19, cex=ce<-0.4)
  cr()
  points(u[1],u[2], col=2)
}
}


# 3d:
if(d3){
library(rgl)
xs <- seq(-1,1, length=20)
x <-as.matrix( expand.grid(x=xs, y=xs, z=xs) )
co <- 1:length(x)
plot3d(x, aspect=FALSE, size=5, col=co)

azis <- seq(0, pi, length=10)[-1]
e <- length(azis):1
while(1)
for(i in e<-rev(e)){
  rgl.clear()
  u <- 0.8*ai2xyz(cbind(0, azis[i]))[1,]#c(0,0,0.8)
  #xx <- xyzrotate(xyzrotate(x, ay=azis[i]), az=azis[i])
  xx <- compress(x, u, T, F)
  points3d(xx, col=co, size=5+i, alpha=0.5)
  lines3d(rbind(0,0,0,u), col=i)
  Sys.sleep(0.1)
}
}