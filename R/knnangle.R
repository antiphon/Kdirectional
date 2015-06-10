#' k-th nearest neighbour angles and distances
#' 
#' @param x Matrix of coordinates, n rows and d columns
#' @param k k-th nearest neighbour
#' @param from From these indices (default all 1,...,nrow(x))
#' @param to To these indices (default all)
#' If 2d, returns the angle in [-pi,pi]. In 3d also the angle from xy-plane in
#' [0,pi]
#' 
#' The values in 2d are (azimuth, distance)
#' The values in 3d are physical coordinate system (azimuth, inclination, distance)
#' 
#' @export

knnangle <- function(x, k=1, from, to){
  x <- as.matrix(x)
  if(missing(from)) from <- 1:nrow(x)
  if(missing(to)) to <- 1:nrow(x)
  d<-c_knnangles(x, k, from, to)
  v<-do.call(cbind, d)
  colnames(v) <- c("ang", if(ncol(v)==2) "dist" else c("ang2", "dist"))
  class(v) <- c("nnangle", is(v))
  attr(v, "k") <- k
  v
}


