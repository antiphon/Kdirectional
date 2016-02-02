#' Test the Kest_anin  in 3D

library(devtools)
load_all(".")

# 
if(!exists("zs")){
  library(rstrauss)
  cat("simulating")
  z  <- zs <- rstrauss(1000, 0.01, 0.05, perfect=T, iter=1e5, bbox=cbind(0:1,0:1,0:1)-.5)
}
bbox <- z$bbox
X0 <- X <- z$x
r <- seq(0, 0.1, length=50)
###############################################
# check the stationary isotropic case
if(0){
  lambda0 <- nrow(z$x)
  out <- Kest_anin(z, epsilon=pi/2, lambda=rep(lambda0, nrow(z$x)))
  library(spatstat)
  pp <- pp3(X[,1],X[,2],X[,3], bbox)
  k0 <- K3est(pp, rmax=max(out$r))
  plot(k0)
  lines(out$r, out$`(1,0,0)`, col=5, lwd=2)
}
### Seems to work.

###############################################
## Check inhomogeneous isotropic case
if(0){
  X <- X0[runif(nrow(X0)) < 0.5+X0[,1],]
  bwv <- seq(0.1, 1, l=20)
  loss <- intensity_bandwidth_profile(X, bwv)
  bw <- bwv[which.min(loss)]
  lambda <- intensity_at_points(X, bw)
  out<- Kest_anin(list(x=X, bbox=bbox), lambda=lambda, r=r)
  oute<- Kest_anin(list(x=X, bbox=bbox), lambda=rep(nrow(X), nrow(X)), r=r)
  plot(out)
  lines(oute[,-2])
}
# works sort of
###############################################
## Inhomogeneous and anisotropic. No real baseline here...
#
## 1. Poisson simulation
if(0){
  vz<-vx <- vy <- NULL
  r <- seq(0, .3, length=50)
  for(i in 1:100){
    X <- matrix(runif(200*3, -.5,.5), ncol=3)
    n <- nrow(X)
    p <- X[,1]
    p <- (p-min(p))/(max(p)-min(p))
    keep <- runif(n) < p
    X <- X[keep,]
    out<- Kest_anin(list(x=X, bbox=bbox), lambda_h=0.7, r=r, epsilon=pi/4)
    vx <- rbind(vx, out$`(1,0,0)`)
    vy <- rbind(vy, out$`(0,1,0)`)
    vz <- rbind(vz, out$`(0,0,1)`)
  }
  the <- out$theo
  
  plot(r,the)
  zx <- apply(vx, 2, quantile, prob=c(.1,.5,.9))
  zy <- apply(vy, 2, quantile, prob=c(.1,.5,.9))
  zz <- apply(vz, 2, quantile, prob=c(.1,.5,.9))
  lines(r,zx[1,], col=1)
  lines(r,zx[3,], col=1)
  lines(r,zx[2,], col=2)
  lines(r,zy[1,], col=1, lty=2)
  lines(r,zy[3,], col=1, lty=2)
  lines(r,zy[2,], col=2, lty=2)
  lines(r,zz[1,], col=1, lty=3)
  lines(r,zz[3,], col=1, lty=3)
  lines(r,zz[2,], col=2, lty=3)
} # ok works quite well!


# 2. Strauss + independent thinning. Should be the same.
if(0) {
  n <- nrow(z$x)
  
  k1 <- Kest_anin(z, epsilon=pi/4, lambda=rep(n,n), r=r)
  
  vz<-vx <- vy <- NULL
  
  for(i in 1:50){
    zv <- z
    zv$x <- zv$x[ runif(n) < .8+zv$x[,1], ]
    out<- Kest_anin(zv, lambda_h=0.66, r=r, epsilon=pi/4)
    vx <- rbind(vx, out$`(1,0,0)`)
    vy <- rbind(vy, out$`(0,1,0)`)
    vz <- rbind(vz, out$`(0,0,1)`)
  }
  the <- out$theo
  par(mfrow=c(2,2))
  zx <- apply(vx, 2, quantile, prob=c(.1,.5,.9))
  zy <- apply(vy, 2, quantile, prob=c(.1,.5,.9))
  zz <- apply(vz, 2, quantile, prob=c(.1,.5,.9))
  p0 <- function() plot(r,the-the, "l", lty=3, ylim=c(-1,1)*5e-4)
  p0()
  lines(r,zx[1,]-the, col=1)
  lines(r,zx[3,]-the, col=1)
  lines(r, k1$`(1,0,0)`-the, col=3)
  p0()
  lines(r,zy[1,]-the, col=1)
  lines(r,zy[3,]-the, col=1)
  lines(r, k1$`(0,1,0)`-the, col=3)
  p0()
  lines(r,zz[1,]-the, col=1)
  lines(r,zz[3,]-the, col=1)
  lines(r, k1$`(0,0,1)`-the, col=3)
  
  # ok.
}

# 3. Strauss + independent thinning + compression. Should differ.
if(1) {
  par(mfrow=c(2,2))
  n <- nrow(z$x)
  zv <- z
  zv$x <- zv$x[ runif(n) < .8+zv$x[,1], ]
  C <- diag(c(1/sqrt(.5), 1/sqrt(.5), .5))
  zv$x <- zv$x%*%C
  zv$bbox <- bbox_affine(zv$bbox, C)
  k1 <- Kest_anin(z, epsilon=eps <- pi/8, lambda_h=lh<-0.7, r=r)
  k2 <- Kest_anin(zv, epsilon=eps, lambda_h=lh, r=k1$r)
  plot(z$x, asp=1)
  plot(k1)
  plot(zv$x, asp=1)
  plot(k2)
  # ok.
}

