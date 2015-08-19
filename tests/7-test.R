#' test directional graph using cutting inside the Kdirectional


library(Kdirectional)
library(rstrauss)
library(icetools)
# set.seed(21)
n <- 1000
r <- strauss_theory(n, 2)*0.7
p <- rstrauss(n=n, range=r, gamma=1e-9, bbox=cbind(0:1, 0:1))
x <- deform.coords(p$x, dy=0.8)

# 
d0 <- nnangle(p$x)
d<-nnangle(x)
par(mfrow=c(2,1))
e0<-hist(d0[,1])
e1<-hist(d[,1])

# chi2 test
chi2 <- function(x){
  E <- sum(x)/length(x)
  x2 <- sum( (x-E)^2 ) / E
  cat("Test for uniformity\nstatistic:", x2, "\ndf:", length(x)-1, "\np:",
  1-pchisq(x2, df=length(x)-1), "\n")
}

chi2(e0$counts)
chi2(e1$counts)
