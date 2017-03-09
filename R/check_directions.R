#' Sector computation directions
#' 
#' @export

check_directions <- function(theta, epsilon, unit, dim=2, antipode=FALSE, n_dir=12) {
  if(missing(theta) & missing(unit)){
    LIM <- if(antipode) 2*pi else pi
    slice <- LIM/n_dir
    if(missing(epsilon)) epsilon <- slice/2
    theta <- list(ang = seq(epsilon, LIM-epsilon, by=slice) )
    if(dim==3) theta <- append(theta, list(ang2 = seq(epsilon, pi-epsilon, by=slice)))
    unit <- theta_2_unit(theta)
    # grid for computation
    thetagrid <- expand.grid(theta)
    unitgrid <- theta_2_unit(thetagrid)
    
  }
  else if(missing(theta) & !missing(unit) ){
    thetagrid <- theta <- unit_2_theta(unit)
    unitgrid <- unit
  }
  else if(!missing(theta) & missing(unit)){
    unit <- theta_2_unit(theta)
    thetagrid <- expand.grid(theta)
    unitgrid <- theta_2_unit(thetagrid)
  }
  else{
    stop("Something went wrong with direction calculations!!?! Shouldn't happend, contact package creator.")
  }
  
  # grid
  
  if(missing(epsilon)) epsilon <- diff(theta[[1]][1:2])/2
  
  grid <- list(theta=thetagrid, unit=unitgrid)
  
  list(theta=theta, unit=unit, epsilon=epsilon, grid=grid)
}

#' Make a direction grid
#' 
#' @export
direction.grid <- function(n_dir=2, d=2, max=pi/2){
  if(d==3) sphere.grid(n_dir)
  else theta_2_unit(list(seq(0, max, length=n_dir)))
}

#' theta (angle) to unit vector, list format input
#' 
#' Convert from angle to unit vector
#' 
#' @param list with 1 or 2 components, giving the angle(s)
#' 
#' @export

theta_2_unit <- function(theta){
  unit <- cbind( cos(theta[1]), sin(theta[2]))
  if( length(theta)==2 ) unit <- cbind(unit*sin(theta[[2]]), cos(theta[[2]]))
  unit
}

#' Angle to unit vector
#' 
#' @param vector of 1 or 2 components (2d or 3d)
#' 
#' @export

angle_2_unit <- function(angle){
  angle <- cbind(angle)
  unit <- cbind(cos(angle[,1]), sin(angle[,1]))
  if( ncol(angle) == 2 ) unit <- cbind(unit*sin(angle[,2]), cos(angle[,2]))
  unit
}

#' unit vector to theta (angle) list format input
#' 
#' @export

unit_2_theta <- function(unit) {
  if(!is.null(dim(unit))) unit <- cbind(unit)
  dim <- ncol(unit)
  unit <- t( apply(unit, 1, function(u) u/sqrt(sum(u^2))  ) )
  theta <- list(ang = atan2( unit[,2], unit[,1] )  )
  if(dim==3) theta <- list(theta, list(ang2 = acos(unit[, 3])))
  theta
}

#' Unit vector to angle
#' 
#' @param unit vector(s as matrix, cols=dimensions)
#' 
#' @export

unit_2_angle <- function(unit){
  if(!is.null(dim(unit))) unit <- cbind(unit)
  dim <- ncol(unit)
  unit <- t( apply(unit, 1, function(u) u/sqrt(sum(u^2))  ) )
  theta <-  atan2( unit[,2], unit[,1] ) 
  if(dim==3) theta <- cbind(theta, ang2 = acos(unit[, 3]))
  theta
}

#' Radians to degree
#' 
#' @export
rad2deg <- function(rad) rad * 180/pi

#' Degree to radian
#' 
#' @export
deg2rad <- function(deg) deg * pi/180
