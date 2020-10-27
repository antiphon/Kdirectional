#' Directional nearest neighbour distribution
#' 
#' Global/local directed nearest neighbour distance distribution, conical wedge defining direction set.
#' 
#' @param x point pattern
#' @param r ranges to consider, increasing vector 
#' @param direction unit vector(s) of directions of cones
#' @param epsilon Cone central half-angle
#' @param type "global" or "local" (default).
#' @param antipodal Use antipodal symmetry (=double cone)? Default: FALSE
#'
#'@export

nn_directional_distance_distribution <- function(x, r, direction, epsilon, type = "local", antipodal = TRUE){
  typ <- pmatch(type, types <- c("global", "local"))
  if(is.na(typ)) stop("'type' should be 'global' or 'local' ")
  x <- check_pp(x)
  bbox <- x$bbox
  dim <- ncol(bbox)
  sidelengths <- apply(bbox, 2, diff)
  # directions
  if(missing(direction)){ 
    u <- diag(c(1), dim)
  }
  else u <- direction
  #
  # make sure unit vectors
  u <- rbind(u)
  u <- t(apply(u, 1, function(ui) ui/c(sqrt(t(ui)%*%ui))))
  #
  # central half-angle
  if(missing(epsilon)){
    epsilon <- pi/4 
  }
  if(abs(epsilon)>pi/2) stop("epsilon should be in range [0, pi/2]")
  #
  # ranges
  if(missing(r)) {
    lambda <- nrow(x$x) / prod(sidelengths)
    rmax <- sqrt(log(1e+05)/(pi * lambda))
    r <- seq(0, rmax, length=50)
  }
  # setup ready.
  if(typ==1) est <- nn_directional_distance_distribution_global(x, u, epsilon, r, antipodal)
  else if(typ == 2) est <- nn_directional_distance_distribution_local(x, u, epsilon, r, antipodal)
  
  # compile result
  dir_names <- apply(u, 1, function(ui) paste0("(", paste0(round(ui,3), collapse=","), ")" ))
  # theoretical
  theo <- NA
  #
  gest <- data.frame(r=r, theo = theo, est)
  names(gest)[] <- c("r", "theo", dir_names)
  rownames(gest) <- NULL
  attr(gest, "epsilon") <- epsilon
  attr(gest, "fname") <- paste0("G_", types[typ])
  #done
  class(gest) <- c("f_cone", is(gest))
  gest
}


#' Directional nearest neighbour distribution, global version
#' 
#' @details 
#' 
#' Called by 'nn_directional_distance_distribution' function.
#' 
#' @export

nn_directional_distance_distribution_global <- function(x, u, epsilon, r, antipodal){
  dim <- ncol(x$x)
  na <- nnangle(x$x)
  ang <- na[,-dim]
  nnd <- na[,dim]
  # Border correction:
  bd <- bbox_distance(x$x, x$bbox)
  bok <- nnd < bd
  d <- nnd[bok]
  # the cone inclusion:
  nnu <- angle_2_unit(ang[bok])
  dev <- apply(u, 1, function(ui ) acos(nnu %*% ui) )
  if(antipodal) dev <- pmin(dev, pi - dev)
  incone <- dev < epsilon
  # counts per range in the sector:
  counts <- apply(incone, 2, function(ic) sapply(r, function(ri) sum(d[ic] <= ri) ))
  t(t(counts) / colSums(incone))
}



#' Directional nearest neighbour distribution, local version
#' 
#' @details 
#' 
#' Called by 'nn_directional_distance_distribution' function.
#' 
#' @export

nn_directional_distance_distribution_local <- function(x, u, epsilon, r, antipodal){
  if(!antipodal) stop("local version only for double cone (antipodal == TRUE)")
  dim <- ncol(x$x)
  bbox <- x$bbox
  bd <- bbox_distance(x$x, bbox)
  bbox0 <- bbox - bbox[1,]
  nr <- length(r)
  r2 <- c(r, Inf)
  e <- NULL
  for(i in 1:nrow(u)){
    na <- nnangle_cone(x$x, unit = u[i,], theta = epsilon, antipodal = antipodal)
    nnd <- na[,dim]
    # Omit NA: If not antipodal, some don't have neighbour in the window
    has_nn <- !is.na(nnd)
    # New border correction:
    bds <- redenbach_border_term(nnd[has_nn], x$bbox, u[1,], epsilon)
    bok <- apply((bd[has_nn] - bds) >= 0, 1, all)
    if(sum(bok))
      w <- apply(sapply(1:dim, function(k) bbox0[2, k] - 2 * bds[bok, k] ), 1, prod)
    else 
      w <- NA
    d <- nnd[has_nn][bok]
    # counts per range:
    counts <- sapply(r2, function(ri) sum( 1/w[d <= ri] ) )
    e <- cbind(e, counts)
  }
  t( t(e[1:nr,])/e[nr+1,]  )
}



#' The erosion of window with sector
#' 
#' @details 
#' 
#' Called by 'nn_directional_distance_distribution' function.
#' 
#' @export

redenbach_border_term <- function(d, bbox, u, epsilon) {
  # shift the box to 0,0,0
  bbox <- bbox - bbox[1,]
  dim <- ncol(bbox)
  # hmm
  norms <- diag(1, dim)
  angle_u_bbox <- acos( u%*%norms )
  beta <- angle_u_bbox - epsilon
  beta[beta<0]<-1
  dfactor <- cos(beta)
  # need to check this
  t(sapply(d, function(di) di*dfactor))
}




#' Directional Nearest neighbour distribution(old)
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



