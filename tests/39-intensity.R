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
v <- intensity_somewhere(X, loc, bw1)
cat("inten\n")
# 
epai <- epa_integral(loc, z$bbox, bw1)
cat("epa\n")
M <- matrix(v, ncol=length(gx))

Mi <- matrix(epai, ncol=length(gx))



par(mfrow=c(2,3))

plot(bwv, err, "l", ylim=c(0,1))

image(gx, gy, Mi, asp=1, zlim=c(0,1), main="epanech integral")

image(gx, gy, M, asp=1, zlim=c(0,2*nrow(X)))
points(X)

print(range(M))
print(c(n=nrow(X),int=sum(M*dx)))

# Compare to Cross validation
library(spatstat)
pp <- ppp(X[,1], X[,2], window = as.owin(c(bbox_make(X))))
t1 <- system.time(bwp <- bw.ppl(pp, srange = range(bwv)) ) # for gaussian

plot(bwp)
image(int2 <- density(pp, sigma = bwp))

print(rbind(t0,t1))

# # use it in an estimate of K
# l1 <- intensity_at_points(X, bw1)
# l2 <- int2[pp]
# k1 <- Kinhom(pp, lambda = l1)
# k2 <- Kinhom(pp, lambda = l2)
# plot(k1$r, k1$trans, type="l")
# lines(k2$r, k2$trans, col=2)
