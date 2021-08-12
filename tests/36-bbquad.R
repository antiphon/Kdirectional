# quad bb
library(devtools)
load_all(".")
library(rgl)

#load("test_strauss_3d.rda")
xc<-matrix(runif(300,-.5,.5) , ncol=3)
bb <- bbq <- bbquad_default(c(-.5,.5))

# rotate
R <- rotationMatrix3(ax=pi/4, az=pi/4)
bbR<-bbox_affine(bb, R)
xcR <- xc%*%t(R)

# compress rotation
v<-0.8
M <- diag(c(1/sqrt(v),1/sqrt(v),v))
V <- M%*%R
xcCR <- xc%*%t(V)
bbz <- bbox_affine(bb, V)
# shift for plotting
sh <- c(3,1,0)
xcCR <- t(t(xcCR)+sh)
bbz <- bbox_affine(bbz, diag(1,3), sh)

# plot
# plot3d(xc, aspect=FALSE)
# plot(bb, alpha=0.1, edges=TRUE)
# 
# points3d(xcCR)
# plot(bbz, alpha=0.1, col=2, edges=TRUE)

# distances: start with the simple bbox
l <- c(-1,1)*.5
bb0 <- cbind(l,l,l)
bd0 <- bbox_distance(xc, bb0)

bd1 <- bbquad_distance(xc, bb)
bdq <- bbquad_distance(xcR, bbR)

cbind(bd0, bd1, bdq)
all.equal(bd0, bd1)  & isTRUE(all.equal(bd0, bdq))

###########################################################

# Fry points with the new bbox thing
#
f0 <- fry_points(pp0<-list(x=xc, bbox=bb0))
f1 <- fry_points(pp1<-list(x=xc, bbox=bb))
f2 <- fry_points(pp2<-list(x=xcR, bbox=bbR))
all.equal(f0$fry_r, f1$fry_r)
all.equal(f0$fry_r, f2$fry_r)
all.equal(f0$fry_units, f1$fry_units)
all.equal(f0$fry_units, f2$fry_units%*%R, check.attributes = F)
## ok!

# fry ellipsoids
ff0 <- fry_ellipsoids(pp0, eps=pi/10)
ff1 <- fry_ellipsoids(pp1, eps=pi/10)
ff2 <- fry_ellipsoids(pp2, eps=pi/10)

plot(ff0$ellipsoids[[3]], aspect=FALSE, col=2)
plot(ff2$ellipsoids[[3]], aspect=FALSE, col=3)


## and then the profiling
co <- seq(0.5,1.1, by=0.1); so <- 1/sqrt(co)
grid <- cbind(so,so,co)
ppz <- list(x=xc%*%M, bbox=bbox_affine(bb, M))
pp0 <- list(x=xc%*%M, bbox=bb0%*%t(M))
pp00 <- list(x=xc, bbox=bb0)
ap<-anisotropy_profile(ppz, grid, r=seq(0,0.2, length=40))
ap0<-anisotropy_profile(pp0, grid, r=seq(0,0.2, length=40))
ap00<-anisotropy_profile(pp00, grid, r=seq(0,0.2, length=40))
plot(ap$profile)
points(ap0$profile,col=2)






