#' test directional graph using cutting

library(Kdirectional)
library(rstrauss)
set.seed(1)
p <- rstrauss(n=10000, range=0.01, gamma=0.001, bbox=cbind(0:1, 0:1))

u <- c(0,1,0)
# g0 <- directed_geom(p$x, u=u, theta=pi/4, r=0.15)
# 
# super <- geom(p$x, r=0.15)
# g1 <- directed_geom_by_cut(p$x, u=u, theta=pi/4, r=0.15, pregraph=super)
# all.equal(g1,g0)

# timing.


t0 <- system.time(g0 <- directed_geom(p$x, u=u, theta=pi/4, r=0.015))

super <- geom(p$x, r=0.02)

t1 <- system.time(g1 <- directed_geom(p$x, u=u, theta=pi/4, r=0.015, pregraph=super))
print(all.equal(g1,g0))
print(rbind(t0,t1))
