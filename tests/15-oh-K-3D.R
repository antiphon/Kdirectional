#' oh-K-3d
#' Ohser-Stoyan

library(devtools)
load_all(".")

library(rstrauss)
x <- rstrauss(gamma=0.01, range=0.07, perfect=T)
x0 <- rstrauss(gamma=1, range=0.07, perfect=T)
est <- osK(x)
est0 <- osK(x0)

k<-Kest_directional(x, u=c(0,0,1), theta = pi/2)
plot(k$r, k$trans, "l", col=3)
lines(est$K_iso)
