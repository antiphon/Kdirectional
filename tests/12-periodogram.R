#' test periodogram, Bartlett

library(devtools)
load_all(".")

if(!exists("nclust")){
library(rstrauss)
library(spatstat)

y <- rstrauss(gamma=1, bbox=cbind(0:1, c(0,1)))
x <- rstrauss(gamma=0.01, range=0.05, bbox=y$bbox, perfect=T)
nclust <- function(x0, y0, radius, n) runifdisc(n, radius, centre=c(x0, y0))
z <- list(x=as.matrix(coords(rNeymanScott(10, 0.1, nclust, radius=0.1, n=10, win = as.owin(c(y$bbox))))), bbox=y$bbox)
}

# Bartlett's version
if(1){
per <- function(x) periodogram(x, p=35, std=F)
est0 <- per(y)
est <- per(x)
est1 <- per(z)

zlim <-range(est$Z, est0$Z, est1$Z, na.rm=T)

par(mfrow=c(1,3))
im <- function(x) with(x, image(x=omega[,1], y=omega[,2], z=Z, asp=1, zlim=zlim))
im(est0)
im(est)
im(est1)
}
# Mugglestone Renshaw 2001
if(0){
per <- function(x) {
  d <- DFT(x, m=20)
  d$f <- d$A^2+d$B^2
  d
}
est0 <- per(y)
est <- per(x)
est1 <- per(z)

zlim <-range(est$f, est0$f, est1$f, na.rm=T)

par(mfrow=c(2,3))
im <- function(x) with(x, image(x=pv, y=qv, z=matrix(f, ncol=length(qv)), asp=1, zlim=zlim))
im(est0)
im(est)
im(est1)
#plo <- function(v) with(v, plot(sqrt(rowSums(omegagrid^2)), v, type="p", cex=.2))#, ylim=zlim))
plo <- function(v) with(v, plot(sqrt(rowSums(grid^2)), f, type="p", cex=.2))#, ylim=zlim))
plo(est0)
plo(est)
plo(est1)
}