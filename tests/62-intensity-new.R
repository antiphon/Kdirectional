# Test intensity estimation, new version
library(devtools)
load_all(".")

########################################
# 2D: basic
if(0){

load("test_strauss.rda")

# thin
xx <- z$x
n <- nrow(xx)
p <- xx[,1]
p <- (p-min(p))/(max(p)-min(p))
keep <- runif(n) < p
X <- xx[keep,]
## Optimise best bw
t0 <- system.time(err <- intensity_bandwidth_profile(X))

# estimate on data
intd <- intensity_at_points(X, bw = err$opt, kernel = err$kernel, loo=F)
intd0 <- intensity_somewhere(X, X, bw = err$opt, kernel = err$kernel)

# estimate on a grid
gs <- seq( -.5, .5, l=50 )
grid <- as.matrix(expand.grid(gs,gs))
intg <- intensity_somewhere(X, grid, bw = err$opt, kernel = err$kernel)

# should be same ball park
print(cbind(d=range(intd), grid=range(intg)))

# plot
par(mfrow=c(2,2))

plot(err)
plot(intd-intd0)

zlim <- range(intg, intd)
cm <- heat.colors(120)
ci <- round( (intd-zlim[1])/diff(zlim) * 120) + 1
image(gs, gs, matrix(intg, nrow=50), asp=1, col = cm, zlim = zlim)
points(X, col = cm[ci], pch=19)

}


################################################
# 3D:
if(0){
  load("test_strauss_3d.rda")
  # thin
  xx <- x
  n <- nrow(xx)
  p <- xx[,1]
  p <- (p-min(p))/(max(p)-min(p))
  keep <- runif(n) < p
  X <- xx[keep,]
  ## Optimise best bw
  t0 <- system.time(err <- intensity_bandwidth_profile(X))
  
  # estimate on data
  intd <- intensity_at_points(X, bw = err$opt, kernel = err$kernel, loo=F)
  intd0 <- intensity_somewhere(X, X, bw = err$opt, kernel = err$kernel)
  
  # estimate on one of the planes
  xs <- seq(-.5,.5, l = 50)
  
  grid <- as.matrix(expand.grid(xs,0,xs))
  
  intg <- intensity_somewhere(X, grid, bw = err$opt)
  # thin slice of points
  eps <- 0.05
  d <- 2
  near0 <- X[,d] < eps & X[,d]> -eps
  x2 <- X[near0,-d]
  
  # plot
  par(mfrow=c(2,2))
  
  plot(err)
  plot(intd-intd0)
  
  zlim <- range(intg, intd)
  cm <- heat.colors(120)
  ci <- round( (intd[near0]-zlim[1])/diff(zlim) * 120) + 1
  image(gs, gs, matrix(intg, nrow=50), asp=1, col = cm, zlim = zlim)
  points(x2, col = cm[ci], pch=19, cex=2)
  points(x2, col = "gray50", pch=".")
  
  if(0){
    cm <- heat.colors(120)
    ci <- round( (intd-zlim[1])/diff(zlim) * 120) + 1
    library(rgl)
    plot3d(X, aspect=F)
    spheres3d(X,  color = cm[ci], radius = 0.02)
  }
}

# inputs 2d
if(1){
  load("test_strauss.rda")
  # thin
  X <- z
  ## Optimise best bw
  t0 <- system.time(err <- intensity_bandwidth_profile(X, kernel = "epa"))
  
  
  
  # estimate on data
  intd <- intensity_at_points(X, bw = err$opt, kernel = err$kernel, loo=F)
  
}