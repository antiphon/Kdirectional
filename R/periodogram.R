#' Bartlett's periodogram
#' @param x point pattern, list of $x and $bbox
#' @param p number of frequencies, -p:p x dim grid
#' @param pgrid the grid of integer spectral frequencies
#' @param std Standardize coordinates?
#' @param ... ignored
#' 
#' Omega=0 value is set to NA, and returned as $zerov in the return object.
#' 
#' @export

periodogram <- function(x, p, pgrid, std=TRUE, ...) {
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  bbox <- x$bbox
  dim <- ncol(bbox)
  sidelengths <- apply(bbox, 2, diff)
  n <- nrow(x$x)
  if(missing(p)) p <- n/2
  #'
  #' The grid of integral values
  if(missing(pgrid)){
    pv <- matrix( rep(-p:p, dim), ncol = dim)
    pgrid <- as.matrix( expand.grid(data.frame( pv )) )
    p <- nrow(pv)
  }
  else{
    pv <- apply(pgrid, 2, unique)
    p <- if(is.list(pv)) length(pv[[1]]) else nrow(pv)
  }
  
  omega <- sapply(1:ncol(pgrid), function(i) unique(pgrid[,i])*2 * pi/sidelengths[i])
  omegagrid <- 2*pi*t(t(pgrid)/sidelengths)
  #'
  #'
  loc <- x$x 
  #' "Standardize" coordinates
  if(std){
    loc <- n * t( t( loc/sidelengths ) )
  }
  #'
  #' Compute:
  CONST <- (1i) * 2 * pi/n
  g <- function(pvec) sum( exp( CONST * loc%*%pvec ) )
  J <- apply(pgrid, 1, g)
  #' The constant is chosen mu A = n
  est <- (Re(J)^2 + Im(J)^2) * 2/n
  #'
  #' Drop (0,0,..) value
  zero <- which( apply(pgrid==0, 1, all) )
  zerov <- est[zero]
  est[ zero ] <- NA
  
  list(v=est, Z=matrix(est, ncol=p ), 
       omega=omega, omegagrid=omegagrid,  
       pv=pv, pgrid=pgrid, 
       std=std, bbox=bbox, zerov=zerov)
}


