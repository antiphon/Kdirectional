#' test profile
library(Kdirectional)
library(icetools)

comp <- 0.7
x <- rstrauss(n=200, gamma=0.01, range=0.05, bb=cbind(c(0,1), c(0,2)) )
y <- list(x=deform.coords(x$x, dy=comp), bbox=deform.coords(x$bb, dy=comp))

prof <- anisotropy_profile(y, i=2, cvec=3:11/10, theta=pi/4, r=seq(0.01, 0.1, length=20))
plot(prof$profile, type="b")
abline(v=comp)
