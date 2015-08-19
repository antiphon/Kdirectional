# compare average ellipsoid keeping data and simulting data
#
library(devtools)
library(sphere)
load_all(".")

if(!exists("pp"))source("test/test-data.R")


ff <- fry_ellipsoids(pp, nvec=1:50, keep_data=FALSE, r_adjust=0.25)
ff2 <- fry_ellipsoids(pp, nvec=1:50, keep_data=TRUE, r_adjust=.25)
ave0 <- mean_ellipsoids(ff$ellipsoids, just_rotation = FALSE)
ave <- mean_ellipsoids(ff$ellipsoids, just_rotation = TRUE)
ave2 <- mean_ellipsoids(ff2$ellipsoids, just_rotation = TRUE, keep_data=TRUE)



plot(ff)
plot(ave0, col=1, lwd=5)
plot(ave, col=2, lwd=5)
plot(ave2, col=3, lwd=5)
points(ave2$data, pch=".", col=2)


a1<-sapply(ff$el, getElement, "rot_angle")
a2<-sapply(list(ave0,ave,ave2), getElement, "rot_angle")
plot(a1,1:length(a1))
points(a2, 25:27, col=1:3, pch=19)
