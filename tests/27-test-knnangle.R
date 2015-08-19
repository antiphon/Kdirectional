#' test knnangle
x<-matrix(runif(50*50), ncol=2)
library(devtools)
load_all(".")

nn <- nnangle(x)
vv <- NULL
kvec <- 2:5*4
for(k in kvec) vv<-cbind(vv, knnangle(x, k)[,2])
D <- as.matrix(dist(x))
Dnn <- t(apply(D,1,sort, decreasing=FALSE))[,kvec+1]

print( all(Dnn==vv) )
