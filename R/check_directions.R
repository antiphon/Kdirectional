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
    #' grid for computation
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
  
  #' grid
  
  if(missing(epsilon)) epsilon <- diff(theta[[1]][1:2])/2
  
  grid <- list(theta=thetagrid, unit=unitgrid)
  
  list(theta=theta, unit=unit, epsilon=epsilon, grid=grid)
}


#' theta (angle) to unit vector
#' 
#' @export

theta_2_unit <- function(theta){
  if(!is.list(theta)) stop("theta should be a list.") 
  unit <- cbind( cos(theta[[1]]), sin(theta[[1]]))
  if( length(theta)==2 ) unit <- cbind(unit*sin(theta[[2]]), cos(theta[[2]]))
  unit
}

#' unit vector to theta (angle)
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

