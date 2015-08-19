#' test directional graph using cutting inside the Kdirectional


library(Kdirectional)
library(rstrauss)
set.seed(1)
p <- rstrauss(n=1000, range=0.01, gamma=0.001, bbox=cbind(0:1, 0:1))

u <- c(0,1,0)
# timing.

r<-seq(0,.3, length=50)
t0 <- system.time(g0 <- Kest_directional(p,u,r=r, theta=pi/4))

super <- geom(p$x, r=0.3)

t1 <- system.time(g1 <- Kest_directional(p,u,r=r, theta=pi/4, pregraph=super))
print(all.equal(g1,g0))
print(rbind(t0,t1))
