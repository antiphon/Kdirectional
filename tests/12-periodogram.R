#' test periodogram, Bartlett

library(devtools)
load_all(".")

library(rstrauss)
library(spatstat)

y <- rstrauss(gamma=1, bbox=cbind(0:1, c(0,1)))
x <- rstrauss(gamma=0.01, range=0.05, bbox=y$bbox, perfect=T)
nclust <- function(x0, y0, radius, n) runifdisc(n, radius, centre=c(x0, y0))
z <- list(x=as.matrix(coords(rNeymanScott(10, 0.1, nclust, radius=0.1, n=10, win = as.owin(c(y$bbox))))), bbox=y$bbox)

per <- function(x) periodogram(x, p=19)
est0 <- per(y)
est <- per(x)
est1 <- per(z)

zlim <-range(est$Z, est0$Z, est1$Z, na.rm=T)

par(mfrow=c(1,3))
im <- function(x) with(x, image(x=omega[,1], y=omega[,2], z=Z, asp=1, zlim=zlim))
im(est0)
im(est)
im(est1)
