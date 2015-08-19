#' Check translation vector computation

library(devtools)
load_all(".")

library(spatstat)

set.seed(1)
if(!exists("p")) p <- rStrauss(100, 0.1, 0.09)

x <- list(x=coords(p), bbox=cbind(0:1, 0:1))


k<-Kest_directional(x, epsilon=pi/2, u=c(0,1))
if(!exists("K"))K <- Kest(p, r=k$r)

#'
w <- translation_weights(x)
e <- pairwise_vectors(x$x)
trans <- NULL
for(r in k$r){
  i <- e[,1] < r
  v <- 2*sum(1/w[i])
  trans <- c(trans, v)
}
trans <- trans/(nrow(x$x))^2
#'
plot(k$r, k$trans, type="l", lwd=9)
lines(K$r, K$trans, type="l", col=2, lwd=3)
lines(k$r, trans, col=5, lwd=2)
