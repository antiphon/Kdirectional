#' Estimate the nearest neighbour distance distribution
#' 
#' Spatstat notation G
#' 
#' @export

Gest <- function(x, r){
  warning("Not properly tested.")
  x <- check_pp(x)
  bbox <- x$bbox
  dim <- ncol(bbox)
  #' compute nnearest neighbours 
  a <- nnangle(x$x) 
  nd <- a[,dim] # last column are the distances
  # border correction
  bd <- bbox_distance(x$x, bbox)
  d <- nd[nd < bd] 
  #' cdf
  if(missing(r)) r <- seq(0, max(d)*1.01, length=100) 
  e<-data.frame(r=r, v=ecdf(d)(r))
  #' theo
  lambda <- nrow(x$x)/bbox_volume(bbox)
  A <- if(dim==3) 4/3*pi*r^3 else pi*r^2
  e$theo <- 1-exp(-lambda * A)
  #' done
  e
}