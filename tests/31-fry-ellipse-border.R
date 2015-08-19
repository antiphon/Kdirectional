#
library(devtools)
library(sphere)
load_all(".")

#source("test/test-data.R")
ff <- fry_ellipsoids(pp, nvec=1:50, )
# f0 <- fry_ellipsoids(pp, nvec=nv<-(1:120)*3,r_adjust=ra<-10, double=FALSE)
# f1 <- fry_ellipsoids(pp, nvec=nv, r_adjust = ra, border=TRUE)
# 
# 
# par(mfrow=c(3,1))
# plot(f0)
# plot(f1)
# 
# s0 <- sapply(f0$ellipsoids, getElement, "semi_axes")
# s1 <- sapply(f1$ellipsoids, getElement, "semi_axes")
# plot(s0[1,]/s0[2,], type="b", ylim=c(0,1.2))
# lines(s1[1,]/s1[2,], type="b", col=3)
# abline(h=1, col=2)
