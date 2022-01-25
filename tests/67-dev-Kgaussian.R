# Dev the gaussian kernel second order anisotropy estimator. Idea is to count fry points with anisotropic weights.


# Two issues: Parametrisation and normalisation.
# 
# Martina's parametrisation:

rvec <- seq(0, 30, l = 31)[-1]
r <- 10 # range of the sphere

# Rotation comes from u
R <- t(sphere::rotationMatrix(0,0, -pi/6)[-3,-3])
u <- cbind(1,0) %*% t(R)
d <- ncol(u)
kappa <- .5 # 1/elongation, b/a. Assume 0 < kappa < 1, otherwise mixes with u.
# assume x is main axis, y,z secondary but equal.
C <- diag( c( 1/kappa, rep(kappa^(1/(d-1)), d-1) ) )
S <- R%*%C
# Covariance of a Gaussian so that the ellipsoid covers ~97.5%.
Sig <- S%*%t(S) * (r/1.96)^2
Sig0 <- C%*%C * (r/1.96)^2
#
#
# Projection method for deterimining kernel using isocurves.
ut <- cbind(-u[,2], u[,1]) # 3D: (-u[,2], 0, u[,1]), arbitrary due to b=c in the ellipsoid. "Oblong"? 
# project space spanned by by u, ut, ut2
A <- rbind(u, ut)





# Check, 2D
if(0) {
  # Circle of radius r:
  a <- seq(0, 2*pi, l = 100)#
  x <- cbind(cos(a), sin(a))#
  plot(x*r, asp=1, xlim = c(-1,1)*2*r, ylim = c(-1,1)*2*r, type="l", col = "goldenrod") # circle
  #
  points(u*r, col=4, pch=15) # direction at r-circle 
  # Project to major axis space:
  z <- c(10, 20)
  v <- u*100
  uu <-cbind(-u[2],u[1]) # 2D 
  vv <- uu * 100
  lines(c(0,v[1]),c(0,v[2]), col=4) ;arrows(0,0,vv[1],vv[2], col=4, lty=2) # only 2D plot ok.
  arrows(0,0,z[1],z[2], col=2)
  text(z[1], z[2], "z")
  # length in direction u, for major axis density:
  d1 <- sum(u*z) # dot-product
  points(d1*u, col=z)
  # length in perpendicular direction, any:
  d2 <- sum(uu*z) # this needs the perpendicular vector
  points(uu*d2, col=2)
  d2 <- sqrt( sum(z^2) - d1^2 )
  points(uu*d2, col=3, pch=20)
  if(0){
    # Simulate gaussian:
    # scale range:
    y<- mvtnorm::rmvnorm(5000, c(0,0), sigma = Sig)
    points(y, col=3, cex=.3)
    zy <- y%*%t(solve(S))
    print( mean(rowMeans(zy^2) < r^2) ) # 95%
    #
  }
  points(x%*%t(S) * r, col=5)
}
# Check, 3D:
if(0) {
  library(rgl)
  # arbitrary perpendicular
  u <- c(1,0,0) %*% t(sphere::rotationMatrix(pi/6, pi/7,0))
  l <- c(-10,10)*.2
  plot3d(rbind(0,u), aspect=FALSE, zlim = l, ylim = l, xlim = l)
  lines3d(rbind(0,u))
  z <- cbind(1,1.2,.6)/2
  points3d(z, col="red")
  planes3d(u, alpha=0.2)
  # arbitrary
  uu <- cbind(-u[2],u[1],0)
  planes3d(uu, col=2, alpha=.1)
  # major axis distance
  d1 <- sum(uu*z)
  # any minor axis distance, due to sphericality perpendicular to u
  d2 <- sqrt(sum(z^2) - d1^2)
  points3d(u*d1, col="blue")
  
  points3d(uu*d2, col="blue")
  # should be on a circle in the plane spanned by normal uu.
}

# Quick eval:
if(0) {
  z <- seq(-2,2, l = 200) * r
  gg <- expand.grid(z,z)
  d1 <- c(u%*%t(gg))
  #uu <- cbind(-u[2],u[1]); d2 <- c(uu%*%t(gg))
  d2 <- sqrt(c(rowSums(gg^2) - d1^2))
  va <- dnorm(d1, 0, sqrt(Sig0[1,1])) * dnorm(d2, 0, sqrt(Sig0[2,2]))
  image(z=matrix(va, nrow = length(z), byrow=F), z,z, asp = 1, col = heat.colors(100))
  v <- u*100;lines(c(0,v[1]),c(0,v[2]), col=4)
  a <- seq(0, 2*pi, l = 100)
  x <- cbind(cos(a), sin(a))
  lines(x*r, lty=2)
  points(x%*%t(S) * r, type="l")
  
}


# Ok.



