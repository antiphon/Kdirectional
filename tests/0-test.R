#' 
library(Kdirectional)
library(spatgraphs)

# spherical point cloud
# s <- seq(0, 1, length=20)
# ab <- expand.grid( s*pi,  s*2*pi); a<-ab[,1];b<-ab[,2]
# xx <- sin(a)*cos(b)
# yy <- sin(a)*sin(b)
# zz <- cos(a)
# x <- 0.2 * cbind(xx, yy, zz)
x <- matrix(runif(1000*3, -1,1), ncol=3)
x <- rbind(c(0,0,0), x)
p <- list(x=x[,1], y=x[,2], z=x[,3], n=nrow(x))

inc <- rep(0, nrow(x)); inc[1] <- 1
inc <- inc > 0

R <- 0.9

tt0 <- system.time(g <- g2 <- spatgraph(p, "geo", par=R, include=inc))

tt<-system.time(e <- directed_geom(x, u=c(0,0,1), theta=pi/10, r=R, from=1))
g2$edges <- e
print(tt0)
print(tt)
library(rgl)
plot3d(x, alpha=0.2)
plot(g, p, alpha=0.2, col="blue")
plot.sg(g2, p, col="red")
