# test 2d bbqud
library(devtools)
load_all(".")

b <- bbquad_default(zl=NA)
b3 <- bbquad_default()

a <- pi/3
A <- diag(c(3,2))
R <- matrix(c(cos(a),sin(a),-sin(a),cos(a)),2)

plot(b, add=F,norm=T)

b2 <- bbox_affine(b, A)
plot(b2, add=F,norm=T)

b3 <- bbox_affine(b2, R)

b4 <- bbox_affine(b3, solve(A))

plot(b3, add=F)
plot(b2, ecol=3)
plot(b4, ecol=4)
# check distance
l <- c(-1,1)*.5#bbox_sideLengths(b2)
b0 <- bbox_affine(cbind(l,l),A)

x <- matrix(runif(20),ncol=2)
points(x, cex=.3, col=3)
points(x3<-x%*%t(R), cex=.3, col=2)
points(x4<-x%*%t(R)%*%solve(A), cex=.3, col=4)
d0<-bbox_distance(x,b0)
d2<-bbox_distance(x,b2)
d3<-bbox_distance(x3,b3)
d4<-bbox_distance(x4,b4)


print(all.equal(d0,d2))
print(all.equal(d0,d3))
print(all.equal(d0,d3))


