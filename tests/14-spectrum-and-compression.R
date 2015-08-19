#' Compression

library(devtools)
load_all(".")

library(rstrauss)
library(spatstat)

comp <- 0.6
M <- diag(c(comp, 1/comp))
Mi <- solve(M)
bbox <- cbind(c(0,comp), c(0, 1/comp))
bbox0 <- cbind(0:1, 0:1)

x <- rstrauss(gamma=0.01, range=0.07, bbox=bbox, perfect=T)
nclust <- function(x0, y0, radius, n) runifdisc(n, radius, centre=c(x0, y0))
z <- list(x=as.matrix(coords(rNeymanScott(10, 0.1, nclust, radius=0.1, n=10, win = as.owin(c(x$bbox))))), bbox=y$bbox)

xc <- list(x=x$x%*%Mi, bbox=bbox0)
zc <- list(x=z$x%*%Mi, bbox=bbox0)

#' periodograms
per <- function(x) periodogram(x, p=19)
est <- per(x)
estc <- per(xc)
est1 <- per(z)
est1c <- per(zc)

zlim <- range(est$Z, estc$Z, est1$Z, est1c$Z, na.rm=T)

par(mfrow=c(2,2))
im <- function(x) with(x, image(x=pv[,1], y=pv[,2], z=Z, asp=1))
im(est)
im(estc)
im(est1)
im(est1c)


#' spectrums

per <- function(x) spectrum(x, p=19)
est <- per(x)
estc <- per(xc)
est1 <- per(z)
est1c <- per(zc)




par(mfrow=c(4,2))
im <- function(x, main="") with(x, {
  plot(Rspectrum, type="b", ylim=c(0,3), main=main);abline(h=1)
  plot(Aspectrum, type="b", ylim=c(0,3));abline(h=1)
  abline(v=90)
  })
im(est, "uncompressed Strauss")
im(estc, "Compressed Strauss")
im(est1, "uncompressed N-S")
im(est1c, "Compressed N-S")

#' Does not work for inhibitive??