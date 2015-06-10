#' Nearest neighbour directional distribution
#' 
#' Estimate the directional distribution of Konig&Schmidt 1992
#' with (their notation) 0-nearest neighbour i.e. nearest neighbour.
#' 
#' @param x point pattern, $x coords $bbox bounding box
#' @param r Ranges
#' @param theta Direction sector central angles, list of ncol(x$x)-components.
#' @param epsilon Direction sector width angles
#' @param n_dir If theta not given, greate a grid with this resolution.
#' @details
#' Directions are angle to positive x-axis (2D) or azimuth, inclination (3D). 
#' Theta should be a list of length 1 (2D) and 2 (3D) defining a grid of these angles.
#' 
#' @export

dDirectional <- function(x, r, theta, epsilon, n_dir=10) {
  x <- check_pp(x)
  bbox <- x$bbox
  dim <- ncol(bbox)
  sidelengths <- bbox_sideLengths(bbox)
  #'
  if(missing(theta)) {
    slice <- 2*pi/n_dir
    if(missing(epsilon)) epsilon <- slice/2
    theta <- list(ang = seq(epsilon, 2*pi-epsilon, by=slice) )
    if(dim==3) theta <- append(theta, list(ang2 = seq(epsilon, pi-epsilon, by=slice)))
  }
  if(missing(epsilon)) epsilon <- diff(theta[1:2])/2
  #'
  #' The isotropic case:
  #' compute the nn-angles
  na <- nnangle(x$x)
  #'
  if(missing(r)) {
    b <- max(na[,dim])
    r <- seq(0, b*1.05, length=100)
  }
  #' border:
  bd <- bbox_distance(x$x, x$bbox)
  from <- which(bd > na[,dim])
  #'
  #'
  nr <- na[from, dim]
  G0 <- sapply(r, function(ri) sum(nr < ri)  )
  #' 
  #'
  counts <- NULL
  #' The directed ones
  anggrid <- expand.grid(theta)
  #' make a direction vector out of an angle
  unitv <- function(a){ 
    if(dim==2) c(cos(a), sin(a)) 
    else c( sin(a[2])*cos(a[1]), sin(a[2])*sin(a[1]), cos(a[2]) )
  }
  
  Gt<-apply(anggrid, 1, function(a){
    naa <- nnangle_cone(x$x, unit = unitv(a), theta = epsilon)
    sapply(r, function(ri) sum(naa[from,dim] < ri, na.rm=T)  )
  })
  if(dim==2){
    colnames(Gt) <- paste0("angle_", anggrid[,1])
  }
  else{
    colnames(Gt)  <-  paste0(paste0("azi_", anggrid[,1]), paste0("_inc_", anggrid[,2]))
  }
  est <- Gt/G0
  est[G0==0] <- 0
  list(est=est, r=r, theta=theta, epsilon=epsilon)
}



