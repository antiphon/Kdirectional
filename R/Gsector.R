#' Directional Nearest neighbour distribution
#' 
#' Estimate the directional nearest neighbour CDF.
#' 
#' @param x point pattern, $x coords $bbox bounding box
#' @param r Ranges
#' @param theta Direction sector central angles, list of ncol(x$x)-components.
#' @param unit Alternative to theta, a matrix of direction vectors.
#' @param epsilon Direction sector width angles
#' @param n_dir If theta not given, greate a grid with this resolution.
#' @details
#' 
#' This is a directed sector version of the G-function in spatstat.
#' 
#' Directions are angle to positive x-axis (2D) or azimuth, inclination (3D). 
#' Theta should be a list of length 1 (2D) and 2 (3D) defining a grid of these angles.
#' 
#' @export

Gsector <- function(x, r, theta, unit, directions, epsilon, n_dir=10) {
  x <- check_pp(x)
  bbox <- x$bbox
  dim <- ncol(bbox)
  sidelengths <- apply(bbox, 2, diff)
  #
  if(!missing(directions) & missing(epsilon)) stop("Need epsilon when directions given.")
  if(missing(directions)){
    dirs <- check_directions(theta, epsilon, unit, dim=dim, antipode=TRUE, n_dir=n_dir)
    epsilon <- dirs$epsilon
  }else{
    dirs <- list()
  }
  #
  # The isotropic case:
  # compute the nn-angles
  na <- nnangle(x$x)
  #
  if(missing(r)) {
    b <- max(na[,dim])
    r <- seq(0, b*1.05, length=100)
  }
  # border:
  bd <- bbox_distance(x$x, x$bbox)
  from <- which(bd > na[,dim])
  #
  # 
  #
  counts <- NULL
  # The directions grid
  dirgrid <- if(missing(directions)) dirs$grid$unit else directions
  # make a direction vector out of an angle
  
  Gt<-apply(dirgrid, 1, function(u){
    naa <- nnangle_cone(x$x, unit = u, theta = epsilon)
    ecdf(naa[,dim])(r)
  })

  est <- Gt
  
  A <- if(dim==2) epsilon * r^2 else pi * r^2 * (1-cos(epsilon))
  
  lambda <- nrow(x$x)/prod(sidelengths)
  
  theo <- 1 - exp( - lambda * A)
  
  res <- list(est=est, r=r, directions=dirgrid, 
              unit=dirs$unit, 
              theta=dirs$theta, 
              epsilon=epsilon, theo=theo, dim=dim)
  #
  class(res) <- c("Gsector", "sector", is(res))
  res
}



