#' etst nnangle_cone
library(devtools)
load_all(".")

x <- matrix(runif(10^2),ncol=2)
nc <- nnangle_cone(x, unit=c(1,1), pi/10)


#' plot to make sure

plot(x, asp=1)

l <- function(ar) ar[2]*c(cos(ar[1]), sin(ar[1]))
for(i in 1:nrow(x)){
  n <- nc[i,]
  xy <- rbind(x[i,], x[i,]+l(n))
  arrows(xy[1,1],xy[1,2],xy[2,1],xy[2,2], code = 2, length = 0.1, col=i)
}

# add the cone
u <- attr(nc, "unit")
a <- attr(nc, "theta")
phi <- atan2(u[2],u[1])

for(i in 1:nrow(x)){
  c0 <- x[i,]
  lines(rbind(c0, c0+l(c(phi-a, nc[i,2])), c0+l(c(phi+a, nc[i,2])), c0), col=i)
}
