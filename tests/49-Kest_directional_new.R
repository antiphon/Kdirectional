# test the new Kest directional
library(devtools)
source("tests/test-data.R")
load_all(".") # note order, rstrauss overwrites bbox_distance !
d3 <- 0

x <- x0 <- pp
R <- rotationMatrix3(az=pi/6)[-3,-3]
x2 <- list(x=coord_affine(x$x, R), bbox=bbox_affine(x$bbox, R))
#x2 <- list(x=x$x, bbox=bbox2bbquad(x$bbox))
#plot(bbox2bbquad(x$bbox), add=F); points(x$x, cex=.4)
#plot(x2$bbox, add=F); points(x2$x, cex=.4)
u <- diag(1,2)
ur <- u%*%t(R)
k <- Kest_directional(x, u=u, epsilon=e<-pi/6, border=bo<-1, cylindrical = cyl <- F)
k2 <- Kest_directional(x2, u=ur, e=e, border=bo, cylindrical = cyl)

print(all.equal(bbox_distance(x$x, x$bbox), bbox_distance(x2$x, x2$bbox)))

par(mfrow=c(2,1))
plot(k, rmax=.3)
plot(k2, rmax=.3)


## time comparison to old
t0 <- system.time(k0<-Kest_directional_old(x, u=c(0,1), e=e))
t1 <- system.time(k1<-Kest_directional(x, u=c(0,1), e=e, r=k0$r))
print(rbind(t0, t1))
print(all.equal(unname(k0$trans), k1$`(0,1)`))
