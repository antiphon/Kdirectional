#' mark correlation anin
#' 

library(devtools)
load_all(".")


###### check by homogeneneous
if(0){
set.seed(14)
x <- matrix(runif(400), ncol=2)
marks <- exp(rnorm(nrow(x), x[,1], 0.1))

#plot(x, cex = marks/2, asp=1)

res <- markcorr_anin(x, marks, lambda = nrow(x))
}

# check the pcf
if(0){
resg <- pcf_anin(x, epsilon = pi/4, lambda = nrow(x), divisor = "r")
plot(resg)
lines(res$pcf_anin$r, res$pcf_anin$`(1,0)`)
}
# ok

# check against  spatstat
if(0){
library(spatstat)
pp <- ppp(x[,1],x[,2], marks=marks)
mc <- markcorr(pp)
res <- markcorr_anin(x, r=mc$r, marks, lambda = nrow(x), epsilon = pi/2, stoyan=0.25)
plot(mc)
lines(res$markcorr_anin$r, res$markcorr_anin$`(1,0)`, col=3)
}


### Inhomogeneous
if(0){
keep <- runif(nrow(x)) < x[,2]
xt <- x[keep,]
mt <- marks[keep]
# wrong
res_b <- markcorr_anin(xt, mt, lambda = nrow(xt) )

lambda_v <- intensity_bandwidth_profile(x, bw<-seq(0.1, 1, l = 10))
rest <- markcorr_anin(xt, mt, lambda_h = bw[which.min(lambda_v)] )

plot(rest)
lines(res$markcorr_anin$r, res$markcorr_anin$`(1,0)`, col=4)
lines(res_b$markcorr_anin$r, res_b$markcorr_anin$`(1,0)`, col=5, lty=2)
}
## I guess it works. the thing is so ambiguous anyways.



###### check by homogeneneous 3d
if(1){
  set.seed(5)
  n <- 600
  x <- cbind(matrix(runif(n*2), ncol=2), runif(n, 0,.5))
  
  marks <- 100*exp(rnorm(nrow(x), x[,1], 0.1))
  
  l <- intensity_at_points(x, 0.8, b = 31)
  
  res <- markcorr_anin(x, marks, lambda = l, epsilon=pi/4)
  # thin to make inhomogeneous
  keep <- runif(nrow(x)) < x[,2]
  xt <- x[keep,]
  mt <- marks[keep]
  lt <- l[keep]
  rest <- markcorr_anin(xt, mt, lambda = lt, epsilon=pi/4)
  par(mfrow=c(2,1))
  plot(res)
  plot(rest)
}


