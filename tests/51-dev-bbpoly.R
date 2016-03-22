# dev the bbpoly object
library(rgl)
library(devtools)
load_all(".")

bb <- cbind(0:1, 0:1)-.5
bb3 <- cbind(bb, 0:1-.5)

#### 2D
# x2 <- bbox2bbpoly(bb)
# x2 <- poly2bbpoly(c(0, 1, 1.5, 1.2, 0.5, 0),c(0, 0, 0.5, .9, 1, 1))
# a <- x <- x2 <- bbpoly_affine(x2, diag(5,2))
# b <- bbox2bbpoly(cbind(c(2,7),c(.1,2)))
# 
# plot(x2, add=F, ylim=c(-2,6))
# plot(b, ecol=4, add=T)
# 
# #tr <- bbpoly_simplices(x2, T)
#for(i in 1:length(tr)) plot(tr[[i]], ecol=i, add=T, lty=i) 
#print(bbpoly_volume(x2))

##### 3D
x <- a <- cube <- bbpoly_affine(bbpoly_default(), diag(2.5,3))
b <- cubeShift <- bbpoly_affine(cube, diag(c(1,1,1)),s=c(0.2,0.3,0.1))

#plot(a, add=F)
#plot(b, add=T, ecol=4)

#x <- bbpoly_intersection_claudia(a, b)

#plot(inter, add=F, ecol=3, size=10, norm=T)
#plot(bbpoly_affine(y, diag(1,3), s=c(3,0,0)), add=T, ecol=3, size=10, norm=T)

#x3 <- bbox2bbpoly(bb3)
#x$vertices[1,1] <- x$vertices[1,1]-.2
#x$vertices[5,1] <- x$vertices[5,1]-.2
#x <- bbpoly_affine(x, diag(1.00, 3))
#print(bbpoly_volume(x))

# s <- bbpoly_simplices(x, T)
# 
# plot(x, ecol=NA, alpha=0)
# s <- lapply(s, bbpoly_affine, A=diag(0.8,3), center=T)
# lapply(s, plot ,alpha=0.1, faces=T, ecol=NA, add=T)


######################
# Intersection
cube<-bbpoly_default()
#cubeShift<-bbpoly_affine(cube, diag(c(1,1,1)),s=c(0.2,0.3,0.1))
#plot(cube, add=F)
#plot(cubeShift, add=T, ecol=3)
# 
#inter<-bbpoly_intersection_claudia(cube, cubeShift)
#plot(inter, add=T, ecol=4, size=10)

#transformed cubes
C<-diag(c(1/sqrt(.5), 1/sqrt(.5), .5))
C<-sphere::rotationMatrix(ax=pi/4 )
C<-sphere::rotationMatrix(ax=pi/4 )%*%diag(c(1/sqrt(.5), 1/sqrt(.5), .5))

bbC<-bbpoly_affine(cube, C)
bbshiftC<-bbpoly_affine(bbC, diag(c(1,1,1)),s=c(0.2,0.3,0.1))
inter<-bbpoly_intersection_claudia(bbC, bbshiftC)

# plot(bbC, add=F, ecol=2)
# plot(bbshiftC, ecol=3)
# plot(inter, ecol=4, face=T, alpha=.1)
# print(sapply(list(bbC, bbshiftC, inter), bbpoly_volume))

inter <- bbpoly_intersection_claudia(bbshiftC, cube)
inter2 <- bbpoly_intersection_claudia(cube, bbshiftC)
#plot(cube, add=F, ecol=2)
#plot(bbC, ecol=3)
plot(inter, ecol=4, lwd=100, size=5, add=F)
plot(inter2)



