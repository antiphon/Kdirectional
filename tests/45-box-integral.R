# test box kernel integral calculation

library(devtools)
load_all(".")

## setup
bbox <- cbind(0:1, 0:1)

ss <- seq(0,1, length=50)
ssy <- seq(0,1, length=50)
x <- as.matrix(expand.grid(ss,  ssy))

## integrate at one point, should be 1, check against fine grid
print(c(box_integral(cbind(0.5,.5), bbox, bw=b<-.1), 
        box_integral(cbind(0.5,.5), bbox, bw=b, n=30)))


t0 <- system.time( v1 <- box_integral(x, bbox, bw<-0.1) )
t1 <- system.time( v2 <- box_integral_grid(x, bbox, bw, n=231) )
# 
print(rbind(t0, t1))
# 
M <- matrix(v1, byrow=F, nrow=length(ss))
M2 <- matrix(v2, byrow=F, nrow=length(ss))
Md <- M-M2
# 
image2 <- function(M,...) image(x=ss, y=ssy, M, asp=1, col=gray.colors(200, 0,1), ...)

par(mfrow=c(2,2))
image2(M,zlim=0:1)
image2(M2,zlim=0:1)
image2(Md)
hist(c(Md), 100)
print(range(Md))

