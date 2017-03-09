#' Nearest neighbour directional distribution (new version)
#' 
#' Estimate the directional distribution of Konig&Schmidt 1992
#' with (their notation) 0-nearest neighbour i.e. nearest neighbour, the direction
#' set being a conical wedge of the unit-sphere.
#' 
#' @param x point pattern, $x coords $bbox bounding box
#' @param direction Unit vector direction
#' @param epsilon Direction sector's opening half-angle vector. 
#' @param r Maximum range for of nn-distances to consider. If missing, use max(nn-distances)
#' @param antipodal Make nn-vectors antipodally symmetric? 
#' 
#' @details
#' Compute the probability that a nearest neighbour of a typical point is in direction of the
#' cone given the distance is less than r. 
#' 
#' Input should be one direction and many epsilons OR equal amount of directions and epsilons.
#' 
#' Not using a double cone.
#' 
#' @export

nn_directional_angle_distribution <- function(x, direction, epsilon, r, antipodal = FALSE) {
  x <- check_pp(x)
  bbox <- x$bbox
  dim <- ncol(bbox)
  sidelengths <- bbox_sideLengths(bbox)
  #
  if(missing(direction)) {
    direction <- c(1, rep(0,dim-1))
  }
  if(missing(epsilon)) epsilon <- seq(0, pi/ifelse(antipodal, 2, 1), l = 12)
  #
  # compute the nn-angles
  na <- nnangle(x$x)
  nnd <- na[,dim]
  nnangles <- na[,-dim, drop =FALSE]
  nnunits <- angle_2_unit(nnangles)
  #
  if(missing(r)) {
    r <- max(nnd)
  }
  if(length(r)!=1) stop("r should be a single positive number.") 
  #
  # border correction:
  bd <- bbox_distance(x$x, x$bbox)
  bdok <- bd > nnd
  #
  # For each direction / epsilon:
  u <- rbind(direction)
  u <- u / sqrt(rowSums(u^2))
  eu <- as.matrix( data.frame(epsilon, u, row.names = NULL) )
  #
  # estimate for one sector:
  inrange <- nnd <= r
  nnu <- nnunits[bdok & inrange, , drop  = FALSE ]
  # the angles between nn-vector and sector directions
  dev_one <- function(s){
    acos(nnu %*% s[-1])
  }
  deviat <- apply(eu, 1, dev_one)
  # antipodal?
  if(antipodal) {
    deviat <- pmin(deviat, pi-deviat)
  }
  # The sums
  counts <- rowSums( t( abs(deviat) ) <= eu[,1] )
  # Normalising constant:
  normaliser <- nrow(nnu)
  # estimate:
  est <- counts / normaliser
  # done.
  out <- data.frame(direction = u, epsilon=epsilon, est = est, row.names = NULL)
}






#' Nearest neighbour directional distribution (old version)
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
  #
  if(missing(theta)) {
    slice <- 2*pi/n_dir
    if(missing(epsilon)) epsilon <- slice/2
    theta <- list(ang = seq(epsilon, 2*pi-epsilon, by=slice) )
    if(dim==3) theta <- append(theta, list(ang2 = seq(epsilon, pi-epsilon, by=slice)))
  }
  if(missing(epsilon)) epsilon <- diff(theta[1:2])/2
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
  nr <- na[from, dim]
  G0 <- sapply(r, function(ri) sum(nr < ri)  )
  #
  #
  counts <- NULL
  # The directed ones
  anggrid <- expand.grid(theta)
  # make a direction vector out of an angle
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



