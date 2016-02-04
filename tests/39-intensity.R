# Test intensity estimation
library(devtools)
load_all(".")

load("test_strauss.rda")

# thin
xx <- z$x
n <- nrow(xx)
p <- xx[,1]
p <- (p-min(p))/(max(p)-min(p))
keep <- runif(n) < p
X <- xx[keep,]
## Optimise best bw
bwv <- seq(0.05, 1.5, length=20)

t0 <- system.time(err <- intensity_bandwidth_profile(X, bwv))

cat("opt\n")
bw1 <- bwv[which.min(err)]
print(bw1)
# draw the best
gy <- gx <- seq(-.5, .5, length=50)
dx <- diff(gx[1:2])^2
loc <- as.matrix( expand.grid(gx, gy) )

v <- intensity_somewhere(X, loc, bw1, b = 121) # grid
v2 <- intensity_somewhere(X, loc, bw1, b = 1) # box
v3 <- intensity_somewhere(X, loc, bw1, b = 0) # correct

cat("inten\n")
# 
epai <- epa_integral(loc, z$bbox, bw1, n=91)
epaib <- epa_integral(loc, z$bbox, bw1, n=1)
epaic <- epa_integral(loc, z$bbox, bw1, n=0)

cat("epa\n")
M <- matrix(v, ncol=length(gx))
M2 <- matrix(v2, ncol=length(gx))
M3 <- matrix(v3, ncol=length(gx))
Mi <- matrix(epai, ncol=length(gx))
Mib <- matrix(epaib, ncol=length(gx))
Mic <- matrix(epaic, ncol=length(gx))

par(mfrow=c(3,3))

plot(bwv, err, "l", ylim=c(0,1))
co <- rainbow(120)
image(gx, gy, Mi, asp=1, zlim=c(0,1), col=co, main="epanech integral grid")
image(gx, gy, Mib, asp=1, zlim=c(0,1), col=co, main="epanech integral box")
image(gx, gy, Mic, asp=1, zlim=c(0,1), col=co, main="epanech integral correct")
image(gx, gy, M, asp=1, zlim=c(0,2*nrow(X)), col=co, main="grid")
points(X)
image(gx, gy, M2, asp=1, zlim=c(0,2*nrow(X)), col=co, main="box")
image(gx, gy, M3, asp=1, zlim=c(0,2*nrow(X)), col=co, main="correct")
image(gx, gy, M-M2, zlim=zl <- c(-1,1)*20, asp=1, main="grid-box", col=co)
image(gx, gy, M-M3, zlim=zl, asp=1, main="grid-correct", col=co)

print(range(M-M3))
print(c(n=nrow(X),intgrid=sum(M*dx), intbox=sum(M2*dx), intcorec=sum(M3*dx), ssbox=mean((M-M2)^2), sscor=mean((M-M3)^2)))
