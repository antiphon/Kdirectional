#' Test the Kest_anin 

library(devtools)
load_all(".")

# thinning
if(!exists("z"))load("test_strauss.rda")

u <- cbind(0:1,1:0)
bbox <- z$bbox

###############################################
# check the stationary isotropic case
if(1){
lambda0 <- nrow(z$x)
out <- Kest_anin(z, u, pi/3, lambda=rep(lambda0, nrow(z$x)))

if(!exists("k0")){
  library(spatstat)
  k0 <- Kest(ppp(z$x[,1], z$x[,2], window = as.owin(c(bbox))), r=out$r, correction = "trans")
}

plot(k0$r, k0$trans, type="l")
lines(k0$r, k0$theo, col=2)
lines(out$r, out$theo, col=3, lty=2)
lines(out$r, out$`(1,0)`, col=4, lty=2)
lines(out$r, out$`(0,1)`, col=5, lty=3)
lines(out$r, out$theo, col=5, lty=3, lwd=5)
}
### Works

###############################################
## Check inhomogeneous isotropic case
if(0){
out<- Kest_anin(list(x=X, bbox=bbox), lambda=lambda, normpower = 2, u=u)

if(!exists("k0i")){
  library(spatstat)
  pp <- ppp(X[,1], X[,2], window = as.owin(c(bbox)))
  k0i <- Kinhom(pp, r=out$r, correction = "trans", lambda = lambda, normpower = 2)
}

plot(k0i$r, k0i$trans, type="l")
lines(k0i$r, k0i$theo, col=2)
lines(out$r, out$theo, col=3, lty=2)
lines(out$r, out$`(1,0)`, col=4, lty=2)
lines(out$r, out$`(0,1)`, col=4, lty=2)
}
# works
###############################################
## Inhomogeneous and anisotropic. No real baseline here...
#
## 1. Poisson simulation
if(0){
  vx <- vy <- NULL
  r <- seq(0, .3, length=50)
  for(i in 1:100){
    X <- matrix(runif(200*2, -.5,.5), ncol=2)
    n <- nrow(X)
    p <- X[,1]
    p <- (p-min(p))/(max(p)-min(p))
    keep <- runif(n) < p
    X <- X[keep,]
    out<- Kest_anin(list(x=X, bbox=bbox), lambda_h=.3, n=21, u=u, r=r, epsilon=pi/4)
    vx <- rbind(vx, out$`(1,0)`)
    vy <- rbind(vy, out$`(0,1)`)
  }
  the <- out$theo
  
  plot(r,the)
  zx <- apply(vx, 2, quantile, prob=c(.1,.5,.9))
  zy <- apply(vy, 2, quantile, prob=c(.1,.5,.9))
  lines(r,zx[1,], col=1)
  lines(r,zx[3,], col=1)
  lines(r,zx[2,], col=2)
  lines(r,zy[1,], col=1, lty=2)
  lines(r,zy[3,], col=1, lty=2)
  lines(r,zy[2,], col=2, lty=2)
}
# 2. Strauss + independent thinning. Should be the same.
if(0) {
  par(mfrow=c(2,2))
  n <- nrow(z$x)
  zv <- z
  zv$x <- zv$x[ runif(n) < .8+zv$x[,1], ]
  k1 <- Kest_anin(z, u, epsilon=pi/9, lambda_h=0.5, r=seq(0,0.1, length=50))
  k2 <- Kest_anin(zv, u, epsilon=pi/9, lambda_h=0.5, r=k1$r)
  plot(z$x, asp=1)
  plot(k1)
  plot(zv$x, asp=1)
  plot(k2)
  # ok.
}

# 3. Strauss + independent thinning + compression. Should differ.
if(0) {
  par(mfrow=c(2,2))
  n <- nrow(z$x)
  zv <- z
  zv$x <- zv$x[ runif(n) < .8+zv$x[,1], ]
  C <- diag(c(1/.5, .5))
  zv$x <- zv$x%*%C
  zv$bbox <- bbox_affine(zv$bbox, C)
  k1 <- Kest_anin(z, u, epsilon=eps <- pi/8, lambda_h=lh<-0.7, r=seq(0,0.1, length=50))
  k2 <- Kest_anin(zv, u, epsilon=eps, lambda_h=lh, r=k1$r)
  plot(z$x, asp=1)
  plot(k1)
  plot(zv$x, asp=1)
  plot(k2)
  # ok.
}

