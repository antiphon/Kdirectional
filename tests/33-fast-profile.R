# tes fast profiling
# OBSOLETE dev version. "fast" became default.
#
#
#' test
library(devtools)
load_all(".")


l<-c(-1,1)*0.5
gam <- 0.6
gam2 <- (1+gam)/2
M <- diag(c(gam, gam2, 1/(gam*gam2)))

# library(rstrauss)
#x <- rstrauss(1200, 0.1, strauss_theory(1000, 3)*0.6, bbox=cbind(l,l,l)%*%t(solve(M))  , perfect=T, iter=1e5, CFTP_T0 = 1e4-1)
# 
#x <- x$x%*%t(M)
# 
#save(file="test_strauss_3d.rda", x)
load("test_strauss_3d.rda")
x <- pp <- list(x=x, bbox=cbind(l,l,l))

cvec <- seq(0.3, 1.2, by=0.2) #grid <- as.matrix(expand.grid(a<-seq(0.7, 1.3, by=.1), b<-seq(0.8, 1.2, by=.1)))
#grid <- as.matrix(expand.grid(a<-seq(0.7, 1.3, by=.05), b<-seq(0.8, 1.2, by=.05)))
cvec2 <- cvec
cvecg <- as.matrix(expand.grid(cvec,cvec2))#cbind(cvec, cvec)
grid <- cbind(a=cvecg, c=1/apply(cvecg,1, prod))

r <- seq(0,0.05, length=20)

t1 <- system.time(   p <- anisotropy_profile_fast(pp, grid, r=r)  )
#t0 <- system.time(  p0 <- anisotropy_profile(pp, cvec=cvec, epsilon=pi/4, r=r)  )



par(mfrow=c(2,1))
with(p, plot(profile[,1], profile$anisotropy))
#with(p0, plot(profile))

#print(cbind(t1,t0))

ani <- p$profile[,4]
ani <- ani/max(ani) * 100
anim <- matrix(ani, ncol=length(cvec))
ce <- exp(1-ani/100)*2
co <- values2colors(-ani, n = length(ani), col = heat.colors)
ox <- order(ani, decreasing = T)
plot(grid[ox,-3], cex=ce, col=co[ox], pch=19, asp=1)
#text(grid[,1], grid[,2], round(ani) , cex=.75)
points(gam,gam2, col="green", cex=2)
# 

# #'
# 
