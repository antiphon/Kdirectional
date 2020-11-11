# Check the gaussian  type-K

library(devtools)
load_all(".")

# thinning
if(!exists("z")){
  load("test_strauss.rda")
  #z$x <- cbind( runif(100), runif(100) )
}

u <- cbind(0:1,1:0)
bbox <- z$bbox
r <- seq(0, .3, l = 20)[-1]

###############################################
# check the stationary isotropic case
if(0){
  lambda0 <- nrow(z$x)
  r <- seq(0, .2, l = 20)[-1]
  kappa <- .5
  out <- Kest_gaussian(z, u, kappa = kappa, lambda = rep(lambda0, nrow(z$x)), r = r)
  #out[,2] <- out[,3]
  par(mfrow=c(2,1))
  plot(out, ylim  =c(0,2))
  # By hand:
  if(!exists("a0")|1){
    f <- fry_points(z, border=FALSE)
    v <- f$fry_r * f$fry_units
    sigs <- r * 0.5
    w <- translation_weights(z)
    l2 <- lambda0^2
    a0 <- colSums(sapply(sigs, function(s) dnorm(v[,1], 0, s/kappa) * dnorm(v[,2], 0, s*kappa) / w )) * 2 / l2
    lines(r, a0, ylim = c(0,2), col=4, lty=2)
  }
  #v <- lambda0^2
  #lines(out$r, out$theo, col=3, lty=2)
  #lines(out$r, out$`(1,0)`/v, col=4, lty=2)
  #lines(out$r, out$`(0,1)`/v, col=5, lty=3)
  #lines(out$r, out$theo, col=5, lty=3, lwd=5)
}
### Works

# Check variance
if(0) {
  r <- seq(0, .2, l = 20)[-1]
  kappa <- 1
  o <- NULL
  for(i in 1:1000) {
    z$x <- cbind(runif(100), runif(100))
    lambda0 <- nrow(z$x)
    o1 <- Kest_gaussian(z, u, kappa = kappa, lambda = rep(lambda0, nrow(z$x)), r = r)[,3]
    o <- cbind(o, o1)
  }
  m <- rowMeans(o)
  plot(r, m, type="l", ylim = c(0,2))
  apply(o, 2, lines, x = r, col=rgb(0,0,0,.1))
  lines(r, m , col = 3)
  plot(r, apply(o, 1, var))
}
  
###############################################
## Inhomogeneous and anisotropic. No real baseline here...
#
## 1. Poisson simulation
if(0){
  vx <- vy <- NULL
  r <- seq(0, .3, length=50)[-1]
  for(i in 1:100){
    X <- matrix(runif(200*2, -.5,.5), ncol=2)
    n <- nrow(X)
    p <- X[,1]
    p <- (p-min(p))/(max(p)-min(p))
    keep <- runif(n) < p
    X <- X[keep,]
    kappa <- .5
    out <- Kest_gaussian(list(x=X, bbox=bbox), u, kappa = kappa, lambda_h = .25, r = r)
    
    vx <- rbind(vx, out$`(1,0)`)
    vy <- rbind(vy, out$`(0,1)`)
  }
  the <- out$theo
  
  plot(r, the, ylim = c(0,2))
  zx <- apply(vx, 2, quantile, prob=c(.1,.5,.9))
  zy <- apply(vy, 2, quantile, prob=c(.1,.5,.9))
  lines(r,zx[1,], col=1)
  lines(r,zx[3,], col=1)
  lines(r,zx[2,], col=2)
  lines(r,zy[1,], col=1, lty=2)
  lines(r,zy[3,], col=1, lty=2)
  lines(r,zy[2,], col=2, lty=2)
  # ok, intensity estimation could be better.
}

# 2. Strauss + independent thinning. Should be the same or very similar at least.
if(0) {
  par(mfrow=c(2,2))
  n <- nrow(z$x)
  zv <- z
  zv$x <- zv$x[ runif(n) < .5+zv$x[,1], ]
  kappa <- .5
  k1 <- Kest_gaussian(z, u, kappa = kappa, lambda_h = .15, r = r)
  k2 <- Kest_gaussian(zv, u, kappa = kappa, lambda_h = .15, r = k1$r)
  plot(z$x, asp=1)
  plot(k1, ylim = c(0,4))
  plot(zv$x, asp=1)
  plot(k2, ylim = c(0,4))
  # ok.
}

# 3. Strauss + independent thinning + compression. Should differ.
if(0) {
  n <- nrow(z$x)
  zv <- z
  zv$x <- zv$x[ runif(n) < 1+zv$x[,1], ]
  C <- diag(c(1/.75, .75))
  zv$x <- zv$x%*%C
  zv$bbox <- bbox_affine(zv$bbox, C)
  k1 <- Kest_gaussian(z, u, kappa = kappa, lambda_h = .15, r = r)
  k2 <- Kest_gaussian(zv, u, kappa = kappa, lambda_h = .15, r = k1$r)
  k3 <- pcf_directional(zv, u,  r = k1$r, stoyan = .3, divisor="d")
  
  par(mfrow=c(3,2))
  plot(z$x, asp=1)
  plot(k1, ylim = y<-c(0,2))
  plot(zv$x, asp=1)
  plot(k2, ylim = y)
  plot(k3, ylim = c(0,2.5))
  # ok.
}

# 4. 3. in 3d
if(0) {
  n <- 200
  if(!exists("z3")) z3 <- rstrauss::rstrauss(n = n, bbox = cbind(0:1,0:1,0:1)-.5, gamma = .1, r = 0.1, iter = 1e4, toroidal=TRUE)
  zv <- z <- z3
  zv$x <- zv$x[ runif(n) < .8+zv$x[,1], ]
  cc<-.6; C <- diag(c(1/cc, sqrt(cc), sqrt(cc)))
  zv$x <- zv$x%*%C
  zv$bbox <- bbox_affine(zv$bbox, C)
  kappa <- .25
  k1 <- Kest_gaussian(z,  kappa = kappa, lambda_h = .15, r = r <- seq(0, .9, l = 30)[-1])
  k2 <- Kest_gaussian(zv, kappa = kappa, lambda_h = .15, r = k1$r)
  par(mfrow=c(2,1))
  plot(k1, ylim = c(0,2))
  plot(k2, ylim = c(0,2))
}


# 5. re-formulate the kernel so that x-axis reach about 95% 
if(0){
  library(mvtnorm)
  kappa <- 0.25
  r <- 10
  L <- r*1.2
  xg <- yg <- seq(-1,1,l=100) * L
  x <- expand.grid(xg, yg)
  qq <- 2
  S <- r^2/qq^2 * diag( c(1, kappa)   )
  v <- dmvnorm(x, sigma = S)
  image( matrix(v, 100), x = xg, y=yg, asp = 1, col = hcl.colors(120))
  symbols(0,0, circles = r, inches=F, add=T)
}

# 6. look at the average of compressed strauss
if(1){
  library(rstrauss)
  ka <- 1.5
  C <- diag(c(1/ka, ka))
  bb <- cbind(0:1, 0:1)-.5
  bb <- t(solve(C)%*%t(bb))
  xl <- lapply(1:50, function(a) rstrauss(gamma = 0.01, n = 100, perfect = TRUE, range = 0.05, bbox = bb)$x)
  xlc <- lapply(xl, function(x) t(C%*%t(x)) )
  kl <- lapply(xlc, function(x) Kest_gaussian(x, kappa = .5, lambda = rep(nrow(x),nrow(x)), r = seq(0, 0.3, l = 30)))
  del <- sapply(kl, function(a) a[,4]-a[,3])
  plot(kl[[1]]$r, rowMeans(del), ylim = c(-1,1)*.3)
  apply(del, 2, lines, x = kl[[1]]$r)
}

