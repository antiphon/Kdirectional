# test angles function
library(rgltools)
library(Kdirectional)
library(icetools)

## compute in pieces
comp <- function(p, m=15, buffer=0.1){
  bb <- p$bbox
  res <- matrix(ncol=4, nrow=p$n)
  
  zcut <- seq(bb[1,3], bb[2,3], length=m)
  
  for(i in 1:(m-1)) {
    #' targets
    zrange <- zcut[i:(i+1)]
    pick <- which( p$x[,3]>=zrange[1] & p$x[,3]<=zrange[2] )
    around <- which( p$x[,3]>=(zrange[1]-buffer) & p$x[,3]<=(zrange[2]+buffer) )
    #' 
    T0 <- Sys.time()
    angles <- nnangle(p$x, from=pick, to=around)[pick,]
    bd <- bbox.distance(bb, p$x[pick,])
    affected <- angles[,3]>bd
    #'
    res[pick,] <- cbind(angles, affected)
    cat("\r", i, "/", m-1, "(", format(Sys.time()-T0), ")    ")
  }
  cat("\n")
  colnames(res) <- c("ang1", "ang2", "dist", "border")
  res
}
#### plot
plotit <- function(nn, ...){
  ll <- ai2ll(nn[,1:2])
  llb <- ll[nn[,4]==0,]
  plot(llb[,2:1], pch=".", asp=1, ...)
}
###

n<-250000
bb <- cbind(c(0,3),c(0,3),c(0,50))

# ## uniform
# x <- (apply(bb, 2, function(a)runif(n, a[1],a[2])))
# x3<-round(x,3)
# p <- pp3d(x3, bbox=bb)
# nn3 <- comp(p)
# 
# x2<-round(x,2)
# p <- pp3d(x2, bbox=bb)
# nn2 <- comp(p)
# 
# x1<-round(x,1)
# p <- pp3d(x1, bbox=bb)
# nn1 <- comp(p)
# 
# par(mfrow=c(3,1))
# plotit(nn1, main="round to 1 digit")
# plotit(nn2, main="round to 2 digit")
# plotit(nn3, main="round to 3 digit")
m <- round(0.1*n)
x <- apply(bb, 2, function(a)runif(n, a[1],a[2]))
r <- sample(1:n, n-m)
x[r,] <- round(x[r,], 2)
pn <- pp3d(x, bbox=bb)
nnn <- comp(pn)
plotit(nnn, main="round to 2 digit: 10% noise")


#' 


# strauss?
# library(rstrauss)
# v <- bbox.vol(bb)
# r <- strauss_theory(n/v)*0.7
# p <- rstrauss(n=n, gamma=0.01, iter=1e5, bbox=bb)
# p <- pp3d(p$x, bbox=bb)
# e<-comp(p)
#plotit(p, e)
# regular honeycomp ?






