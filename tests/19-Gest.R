# test Gest
library(spatstat)
library(devtools)
load_all(".")



pp <- rpoispp(100)

x <- list(x=coords(pp), bbox=cbind(0:1, 0:1))

g1 <- Kdirectional::Gest(x)
g2 <- spatstat::Gest(pp, r=g1$r)

plot(g2)
lines(g1$r, g1$v, col=4, lwd=4, lty=3)
lines(g1$r, g1$theo, col=1, lwd=3, lty=2)
