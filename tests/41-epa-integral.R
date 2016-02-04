# test epanechnicov kernel integral calculation
#
#

# NOTE: This actually calls for the box kernel integral when not using grid
################
library(devtools)
load_all(".")

if(!exists("z")) load("test_strauss.rda")


bbox <- cbind(0:1, 0:1)
bbox3 <- cbind(bbox, 0:1)

#x <- cbind(c(0.5,0.05), c(0.5, 0.5))
#x3 <- cbind(x, c(0.5,.5))

ss <- seq(0,1, length=50)
ssy <- seq(0,1, length=50)
x <- as.matrix(expand.grid(ss,  ssy))

#bbox <- cbind(range(ss),range(ssy))

t0 <- system.time( v <- box_integral(x, bbox, bw<-0.5) )
t1 <- system.time( v1 <- epa_integral_grid(x, bbox, bw, 51) )
t2 <- system.time( v2 <- epa_integral_biased(x, bbox, bw) )
t3 <- system.time( v3 <- epa_integral_2d(x, bbox, bw))
print(rbind(t0,t1,t2,t3))

##plot(v[,1:2], asp=1, col=gray(v[,3]/max(v[,3])), pch=19, cex=5)

image2 <- function(M,...) image(x=ss, y=ssy, M, asp=1, col=rainbow(200,start=0, end=1), ...)
#points(x, pch=".")

M <- matrix(v, byrow=F, nrow=length(ss))
M1 <- matrix(v1, byrow=F, nrow=length(ss))
M2 <- matrix(v2, byrow=F, nrow=length(ss))
M3 <- matrix(v3, byrow=F, nrow=length(ss))

par(mfrow=c(3,3))
image2(M,zlim=0:1, main="box approx")
image2(M1,zlim=0:1, main="grid")
image2(M2,zlim=0:1, main="biased")
image2(M3,zlim=0:1, main="correct 2d")
#points(x, pch=".")

D1 <- M1-M
D2 <- M1-M2
D4 <- M1-M3
image2(D1, zlim=zl<-c(-1,1)*2e-1, main="grid-box")
image2(D2, zlim=zl, main="grid-biasd")
image2(D4, zlim=zl, main="grid-correct")
SS <- mean(D4^2)
print(c(ss_gridbox=mean(D1^2),
ss_gridbia=mean(D2^2), ss_bb=mean(D3^2), ss_gridcorre=SS)/SS)
