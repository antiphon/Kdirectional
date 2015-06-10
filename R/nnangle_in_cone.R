#' nearest neighbour angles and distances given that angles are in a cone
#' 
#' @param x matrix of coordinates
#' @param unit direction of the cone
#' @param theta angle or width of the cone
#' @param from indices from which to compute
#' @param to indices in which look for nearest neighbour
#' 
#' The values in 2d are (azimuth, distance).
#' The values in 3d are physical coordinate system (azimuth, inclination, distance)
#' 
#' The angle is inside the cone. NA if no neighbour is found.
#' 
#' @export

nnangle_cone <- function(x, unit, theta, from, to){
  x <- as.matrix(x)
  if(missing(from)) from <- 1:nrow(x)
  if(missing(to)) to <- 1:nrow(x)
  #' make sure unit
  unit <- unit/sqrt(sum(unit^2))
  d<-c_angles_in_a_cone(x, unit, theta, from, to)
  v<-do.call(cbind, d)
  
  #' Drop not found
  dim <- ncol(x)
  M<- diff(range(x))*2
  v[v[,dim]>M, ]  <- NA
  #
  colnames(v) <- c("ang", if(ncol(v)==2) "dist" else c("ang2", "dist"))
  class(v) <- c("nnangle", is(v))
  attr(v, "unit") <- unit
  attr(v, "theta") <- theta
  v
}
