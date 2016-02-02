# Test intensity estimation in 3D
library(devtools)
load_all(".")

if(!exists("z"))load("test_strauss_3d.rda")

# thin
xx <- x
n <- nrow(xx)
p <- xx[,1]
p <- (p-min(p))/(max(p)-min(p))
keep <- runif(n) < p
X <- xx[keep,]
## Optimise best bw
bwv <- seq(0.05, 4, length=20)

t0 <- system.time(err <- intensity_bandwidth_profile(X, bwv))

cat("opt\n")
bw1 <- bwv[which.min(err)]
print(bw1)
# draw the best
plot(bwv, err, "l")

# intensity
ss <- seq(0,1, length=25)
g <- as.matrix(expand.grid(ss,ss,ss))

t1 <- system.time(int <- intensity_somewhere(X, g, bw1))

print(rbind(t0,t1))

#
v <- int
library(rgl)
library(scales)
co <- gray(rev(i<-rescale(v)))

inc <-  abs(g[,3]-.5) < 0.1 

plot3d(g[inc,], col=co[inc], aspect=F, alpha=rev(i[inc]))
