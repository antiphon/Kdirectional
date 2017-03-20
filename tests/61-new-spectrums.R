#' test spectrums, new versions

library(devtools)
load_all(".")

library(rstrauss)
library(spatstat)

if(!exists("z")){
gam <- 0.6

y <- rstrauss(gamma=1, bbox=cbind(c(0,2/gam), c(0,gam*2)))
x <- rstrauss(gamma=0.01, range=0.07, bbox=y$bbox, perfect=T)
nclust <- function(x0, y0, radius, n) runifdisc(n, radius, centre=c(x0, y0))
z <- list(x=as.matrix(coords(rNeymanScott(10, 0.1, nclust, radius=0.1, n=10, win = as.owin(c(y$bbox))))), bbox=y$bbox)


# anisotropy
C <- diag(c(gam, 1/gam))
compr <- function(x) list(x=coord_affine(x$x, C), bbox=bbox_affine(x$bbox, C))
y <- compr(y)
x <- compr(x)
z <- compr(z)
}

per <- function(x) bartletts_spectral_density(x, 
                                              k = seq(-36,36,l=100),#seq(-3000, 3000, l=100),#c(-16:16), 
                                              std=T)
est0 <- per(y)
est <- per(x)
est1 <- per(z)




par(mfrow=c(3,1))
ima <- function(x,...) with(x, {
 plot(im(t(log(sdf_matrix)), xcol = stops[[1]], yrow = stops[[2]]),   ...)
})

ima(est0, main="Poisson")
ima(est, main ="strauss")
ima(est1, main="matern cluster")


#' Does not work for inhibitive??