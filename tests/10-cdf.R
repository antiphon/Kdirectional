#' cdf
library(icetools)


s <- function(n)matrix(runif(n*2), ncol=2)
tr <- NULL
bb <- cbind(0:1, 0:1)
for(i in 1:200){
  x <- s(100)
  ang <- nnangle(x)[,1]
  ang <- ang[ which( bbox.distance(bb, x)  > 0.1  ) ]
  cdf <- angles.cdf(ang, n = n <- 15)$cdf
  tr <- rbind(tr, cdf)
}
cdf <- angles.cdf(ang, n = n)
v0 <- cdf$ang/(2*pi)
plot(cdf$ang, cdf$cdf-v0, ylim=c(-1,1)*.2, "l")
apply(tr, 1, function(y) lines(y-v0, x=cdf$ang, col=rgb(.2,.2,.2,.2)))
lines(cdf$ang, cdf$CI5-v0, col=3)
lines(cdf$ang, cdf$CI95-v0, col=3)


#' type 1 error
p <- sum( apply(tr, 1, function(y)  sum((y > cdf$CI95 | y < cdf$CI5) ) >0 ) ) / nrow(tr)
print(c(type_I_error=p))
