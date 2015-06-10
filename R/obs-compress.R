# #' Volume preserving compression function
# #' 
# #' @param x matrix of coordinates, n rows d cols
# #' @param u vector of compression, see Details.
# #' @param vol_preserve Preserve volume?
# #' @param center Use mass center as the rotation reference point. Otherwise origo.
# #' @details
# #'
# #' If compression vector length is 1, nothing happens.
# #' 
# #' If compression vector length is <1, compression vector denotes the minor axis of elliptical transformation.
# #' If >1, the major axis. 
# #' 
# #' We first rotate the pattern, then compress, then back-rotate. 
# #' This transformation will take place with respect to origo, unless center=T in which case we use the 
# #' mass center as rotation point.
# #'
# #' @import sphere
# #' @export
# compress <- function(x, u, vol_preserve=TRUE, center=FALSE){
#   #' rotate, compress, back-rotate
#   
#   dim <- ncol(x)
#   
#   x0 <- if(center) x0 <- apply(x, 2, mean) else rep(0, dim)
#   
#   x <- t(t(x)-x0)
#   
#   azi <- atan2(u[2], u[1])
#   z <- sqrt(sum(u^2))
#   #' rotation required 
#   azi <- if(azi<0) 2*pi+azi else azi
#   
#   if(z==1) return(x)
#   #' 2D
#   if(dim == 2){  
#     #' we use y-axis as reference
#     R <- matrix(c(cos(azi), sin(azi), -sin(azi), cos(azi)), ncol=2)
#     #' compress
#     M <- if(vol_preserve) diag(c(1/z, z)) else diag(c(1, z))
#   } # 3D
#   else{
#     #' we use z-axis as reference
#     inc <- acos(u[3]/z)  
#     #' rotate
#     Ry <- xyzrotationMatrix(ay=inc)
#     Rz <- xyzrotationMatrix(az=azi)
#     R <- Rz%*%Ry
#     #' compress
#     zs <- if(vol_preserve)  1/sqrt(z) else 1
#     M <- diag(c(zs, zs, z))
#   }
#   x1 <- x%*%R
#   x2 <- x1 %*% M
#   #' back-rotate
#   x3 <- x2 %*% solve(R)
#   
#   #' just some check
# #   par(mfrow=c(4,1))
# #   L <- max(x)*c(-1,1)*2
# #   plot(x, asp=1, xlim=L, ylim=L, col=co<-1:nrow(x))
# #   plot(x1, asp=1, xlim=L, ylim=L, col=co)
# #   plot(x2, asp=1, xlim=L, ylim=L, col=co)
# #   plot(x3, asp=1, xlim=L, ylim=L, col=co)
# #   
#   t(t(x3)+x0)
# }

# Dont use this, not generic enough
