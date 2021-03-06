# negative K est
library(rstrauss)
library(Kdirectional)
z  <- rstrauss(1000, 0.01, 0.05, perfect=T, iter=1e5, bbox=cbind(0:1,0:1,0:1)-.5)
G<-sphere::rotationMatrix(ax=pi/2 )
zvflip <- list(x =z$x%*%t(G), bbox = bbox_affine(z$bbox, G))
k1<-Kest_directional(zvflip, epsilon=pi/4, u=c(0,0,1),  r=seq(0,0.15,0.005))
plot(k1$r, k1$trans, type="l", lwd=2, xlim=c(0,0.15), ylim=c(0,0.004))
k2 <- Kest_anin(zvflip, epsilon=eps <- pi/4, lambda_h=lh<-0.7, b=31, r=seq(0,0.15,0.005), border=1)
plot(k2)
