#' test spectrums, G&R

library(devtools)
load_all(".")

library(rstrauss)
library(spatstat)

y <- rstrauss(gamma=1, bbox=cbind(0:1, c(0,2)))
x <- rstrauss(gamma=0.01, range=0.07, bbox=y$bbox, perfect=T)
nclust <- function(x0, y0, radius, n) runifdisc(n, radius, centre=c(x0, y0))
z <- list(x=as.matrix(coords(rNeymanScott(10, 0.1, nclust, radius=0.1, n=10, win = as.owin(c(y$bbox))))), bbox=y$bbox)


# anisotropy
C <- diag(c(0.8, 1/0.8))
comp <- function(x) list(x=coord_affine(x$x, C), bbox=bbox_affine(x$bbox, C))

y <- comp(y)
x <- comp(x)
z <- comp(z)


per <- function(x) spectrum(x)
est0 <- per(y)
est <- per(x)
est1 <- per(z)




par(mfrow=c(3,2))
im <- function(x) with(x, {
  plot(Rspectrum, type="b", ylim=c(0,3));abline(h=1)
  plot(Aspectrum, type="b", ylim=c(0,3));abline(h=1)
  })

im(est0)
im(est)
im(est1)


#' Does not work for inhibitive??