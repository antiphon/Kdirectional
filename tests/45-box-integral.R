# test box kernel integral calculation

library(devtools)
load_all(".")

if(!exists("z")) load("test_strauss.rda")


bbox <- cbind(0:1, 0:1)

ss <- seq(0,1, length=50)
ssy <- seq(0,1, length=50)
x <- as.matrix(expand.grid(ss,  ssy))

print(c(box_integral(cbind(0.5,.5), bbox, bw=b<-.1),
      box_integral(cbind(0.5,.5,.5), cbind(bbox,0:1), bw=b)))


 t0 <- system.time( v1 <- box_integral(x, bbox, bw<-0.5) )
 t1 <- system.time( v2 <- box_integral_grid(x, bbox, bw, n=31) )
# 
# print(rbind(t0, t1))
# 
 par(mfrow=c(1,2))
 M <- matrix(v1, byrow=F, nrow=length(ss))
 M2 <- matrix(v2, byrow=F, nrow=length(ss))
# 
 image2 <- function(M,...) image(x=ss, y=ssy, M, asp=1, col=gray.colors(200, 0,1), ...)
 image2(M,zlim=0:1)
 image2(M2,zlim=0:1)
# #points(x, pch=".")
# 
