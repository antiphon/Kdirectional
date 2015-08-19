#' test anisitropy.abs
library(Kdirectional)
library(icetools)


x <- rstrauss(n=100, gamma=0.01, range=0.06)

estimates <- Kest_along_axis(x, theta=pi/4)

s1 <- anisotropy.abs(estimates=estimates)
s2 <- anisotropy.abs(x, theta=pi/4)

