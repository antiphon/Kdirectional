# test epanechnicov kernel integral calculation

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

t0 <- system.time( v <- epa_integral(x, bbox, bw<-0.9) )
t1 <- system.time( v1 <- epa_integral_grid(x, bbox, bw, 11) )
t2 <- system.time( v2 <- epa_integral_biased(x, bbox, bw) )

print(rbind(t0,t1,t2))

##plot(v[,1:2], asp=1, col=gray(v[,3]/max(v[,3])), pch=19, cex=5)

par(mfrow=c(2,2))
image2 <- function(M,...) image(x=ss, y=ssy, M, asp=1, col=gray.colors(200,0,1), ...)
#points(x, pch=".")

M <- matrix(v, byrow=F, nrow=length(ss))
M1 <- matrix(v1, byrow=F, nrow=length(ss))
M2 <- matrix(v2, byrow=F, nrow=length(ss))
image2(M,zlim=0:1, main="box approx")
image2(M1,zlim=0:1, main="grid")
image2(M2,zlim=0:1, main="biased")
#points(x, pch=".")

D <- M1-M
image2(D, zlim=c(-1,1), main="grid-box")
print(range(D))
