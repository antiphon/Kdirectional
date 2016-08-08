# Translation weights in a bbquad window

library(devtools)
load_all(".")

d <- 2

x <- matrix(runif(d*100)-.5, ncol=d)
b <- bbquad_default(zl=if(d==2) NA else c(-.5,.5))
x0 <- list(x=x, bbox=b)
# distance to border
bd <- bbox_distance(x, b)

###############################
### Helpers
affine_list <- function(x,A) list(x=x$x%*%t(A), bbox=bbox_affine(x$bbox, A))
plot1 <- if(d==2) function(x,...){
  plot(x$bbox, add=F)
  points(x$x)
} else function(x,...){
  plot(x$bbox, add=F, alpha=.1, edges=TRUE)
  points3d(x$x)
}
###############################
# after rotation
par(mfrow=c(1,2))
R <- sphere::rotationMatrix(az=pi/9)[1:d,1:d]
x1 <- affine_list(x0, R)
b1 <- x1$bbox
plot1(x1)
bd1 <- bbox_distance(x1$x, x1$bbox)
print(all.equal(bd,bd1))



C1 <- diag(c(1,1))
b2 <- bbox_affine(x1$bbox, C1)
bbox_volume(b2)
(sum(bbox_sideLengths(b2))^2/(4*pi))
