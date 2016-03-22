# test bbquad triangulation

library(devtools)
load_all(".")




# 2D
if(0){
b <- bbquad_default(xl=0:1, yl=c(1,5), zl=NA)
x <- b

R <- sphere::rotationMatrix(az=pi/10)[-3,-3]
C <- diag(c(5,1.2))
x <- bbox_affine(b, A<-C%*%R)

s <- bbquad_simplices(x)
print(all.equal(sum(simplex_volume(s)),x$volume))

add_s <- function(x){
  for(k in 1:length(x)){
    v <- x[[k]]
    for(i in 1:3)
      lines(rbind(v,v[1,])[0:1+i,], col=(k-1)*3+i, lty=k)
  }
}
plot(x, add=F, ecol="gray50")
add_s(s)
## 2d ok
}
## 3d
if(1){
library(rgl)
b <- bbquad_default()
R <- sphere::rotationMatrix(az=0.4, ay=.1)
C <- diag(c(1,1,1))
V <- matrix(c(1,2,1,2.3,1.1,0.1,1,1,2), 3)
x <- bbox_affine(b, A<-V%*%C%*%R)
# mess about with the vertices, will not change x$volume!
#x$vertices <- x$vertices+runif(8*3)*sign(x$vertices)
#
s <- bbquad_simplices(x)

add_s <- function(x){
  for(k in 1:length(x)){
    v <- x[[k]]
    for(i in 1:4)
      lines3d(rbind(v,v[1,])[0:1+i,], col=(k-1)*4+i, lwd=k)
  }
}

plot(x, add=F, alpha=0, edges=F)
add_s(s)
print(all.equal(sum(simplex_volume(s)), x$volume))
print(all.equal(sum(simplex_volume(s)), bbquad_volume(x)))
}
