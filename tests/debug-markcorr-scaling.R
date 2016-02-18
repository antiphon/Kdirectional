# debug: marcorr skaling
library(devtools)
load_all(".")
library(rstrauss)
bbox <- rbind(c(0,0,0),c(1,5,5))
#bbox <- rbind(c(0,0),c(1,4))
#set.seed(12)
x <- rstrauss(145, gamma = .1, range =.1, bbox = bbox, perfect=T, iter=1e5)
# attach z-coordinates as marks
m <- x$x[,3]
#m <- rnorm(nrow(x$x))
mcor3 <- markcorr_anisotropic(x, marks=m, adjust=c(1,1))
# plot shows the scaling issue for x and y directions
d <- ncol(bbox)
par(mfrow=c(2,d))
for(i in 1:d) plot(mcor3$r, mcor3$mcor[,i], type="l")
for(i in 1:d) plot(mcor3$r, mcor3$pcf[,i], type="l", ylim=c(0,2))

