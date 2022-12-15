#' Project a Point To Plane, 2D Outcome
#' 
#' Plane defined using a point and a normal
#' 
#' @param xyz points to be projected, xyz-matrix
#' @param normal normal of the plane
#' 
#' We will simply rotate the points xyz and drop z-coordinate.
#' Rotation order: 1) Rz 2) Ry
#' 
#' This is not unique w.r.t. rotation around the normal, be aware.
#' 
#' @export
project_to_plane <- function(xyz, normal){
  normal <- normal / sqrt(sum(normal^2))
  #
  a <- atan2(normal[2], normal[1])
  i <- acos(normal[3])
  # all we need is to straighten the z coordinate:
  Az <- matrix(c(cos(a), -sin(a), 0, sin(a), cos(a),0,0,0,1), 3)
  Ay <- matrix(c(cos(i), 0, sin(i), 0,1,0, -sin(i),0,cos(i)), 3)
  Ay%*%Az%*%normal
  #
  xy <-t( Ay%*%Az%*%t(xyz) )[,1:2]
  attr(xy,"inc") <- i
  attr(xy,"azi") <- a
  xy
}

#' Intersection of a 3d Ellipsoid and a Plane
#' 
#' @param x ellipsoid, with $semi_axes and $rot components
#' @param n normal of the plane
#' @param q a point defining the plane, inside the ellipsoid
#' @param r optional, perpendicular to n
#' 
#' @details 
#' Note that there is no rotation as the result is in an alternative
#' basis. The basis is stored 
#' in the returned ellipse's component $basis3d. Third column is n.
#' 
#' @references 
#' P. P. Klein,"On the Ellipsoid and Plane Intersection Equation," Appl. Math., vol. 3, no. November, pp. 1634-1640, 2012.
#' @export

intersect_ellipsoid_plane <- function(x, n, q, r) {
  if(!is(x, "ellipsoid")) stop("x not an ellipsoid object")
  R <- x$rot
  Ri <- t(R)
  D <- diag(1/x$semi_axes)
  # helper
  Df <- function(a,b) t(D%*%a)%*%(D%*%b)
  dot <- function(a,b) c( t(a)%*%b )
  # orthongonal vectors.
  # random unit to start with
  if(missing(r)) {
    u <- runifsphere(1)
    r <- cross(n, u)
  }
  r <- r/sqrt(sum(r^2))
  s <- cross(n, r)
  s <- s/sqrt(sum(s^2))
  # reorient and shift
  qr <- c( (q-x$center)%*%R )
  nr <- c(n%*%R)
  rr <- c(r%*%R)
  sr <- c(s%*%R)
  # recompute the orthogonals to fullfill eq 7
  Drs <-Df(rr,sr)
  Drr <-Df(rr,rr)
  Dss <-Df(sr,sr)
  om <- 0.5 * atan2(2*Drs, Drr-Dss)
  rrr <- rr
  rr <- cos(om)*rr + sin(om)*sr
  sr <- -sin(om)*rrr + cos(om)*sr
  # transform back the originals
  r <- c(rr%*%Ri)
  s <- c(sr%*%Ri)
  # 
  Dqr <- Df(qr,rr)
  Dqs <- Df(qr,sr)
  Drr <- Df(rr,rr)
  Dss <- Df(sr,sr)
  d <- Df(qr,qr) - Dqr^2/Drr-Dqs^2/Dss
  if(1-d<0) d <-1
  # semi axes
  A <- sqrt(  (1-d)/Drr  )
  B <- sqrt(  (1-d)/Dss  )
  # center
  c1 <- -Dqr/Drr
  c2 <- -Dqs/Dss
  # need change of basis?
  # no clear rotation as basis is not fixed
  Rxy <- NA
  #
  out <- as_ellipsoid(semi_axes=c(A,B), 
                      R=diag(1,2),
                      center = c(c1,c2))
  # 3d: center
  basisr <- cbind(rr,sr,nr)
  basis <- cbind(r,s,n)
  #browser()
  mr <- qr + c( basisr%*%c(c1, c2, 0) ) #qr + c1*rr + c1*sr
  m <- mr %*% Ri + x$center
  out$center3d <- m
  # basis, the plane basis is r,s
  out$basis3d <- cbind(r,s, n)
  out
}



#' Compute the areas of polygons per 2D or 3D triangulation
#' 
triangulation_areas <- function(tri) {
  if(is(tri, "mesh3d")){
    apply(tri$it, 2, function(ijk){
      vert <- tri$vb[,ijk]
      vert <- t(vert[1:3,])/vert[4,]
      pl <- three_points_to_plane(vert)
      px <- project_to_plane(vert, pl[1:3])
      area_of_2d_polygon(px)
    })
  }
  else{ # 2d triangulation, no need to project
    apply(tri$it, 2, function(ijk){
      vert <- tri$vb[,ijk]
      px <- if(nrow(vert)==3) t(vert[1:2,])/vert[3,] else t(vert)
      area_of_2d_polygon(px)
    })
  }  
}

#' values to colors
#' 
#' @param v the values
#' @param n number of colors
#' @param zlim limits
#' @param color function, e.g. heat.colors, gray.colors
#' 
values2colors <- function(v, n=100, zlim, col=heat.colors, na.col="gray50", ...){
  zlim0 <- range(v, na.rm = TRUE)
  if(missing(zlim)) zlim <- zlim0
  pr <- v
  na <- which(is.na(pr))
  pr[na] <- mean(pr, na.rm=TRUE)
  pr <- pmax(zlim[1], pmin(zlim[2], pr))
  pp <- (pr-zlim[1])/diff(zlim)
  e <- floor(pp*(n-1)) + 1
  hc <-col(n)
  out <- hc[e]
  out[na] <- na.col
  out
}


#' Compute the area of 2D polygon given by ordered border points
#' 
#' So called shoelace method.
#' 
#' @param px polygon coordinates, open ended counter clockwise
#' 
#' @export
area_of_2d_polygon <- function(px){
  #px <- rbind(px,px[nrow(px),])
  px <- rbind(px,px[1,])
  n <- nrow(px)
  0.5 * abs( sum(px[-n,1]*px[-1,2]) - sum(px[-1,1]*px[-n,2]) )
}


