#' Bartlett's periodogram
#' 
#' @param x point pattern, list of $x and $bbox
#' @param p number of frequencies, (-p:p)^2 grid
#' @param pgrid the grid of integer spectral frequencies, two column matrix.
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
  #
  # The grid of integral values
  if(missing(pgrid)){
    pv <- matrix( rep(-p:p, dim), ncol = dim)
    pgrid <- as.matrix( expand.grid(data.frame( pv )) )
    p <- nrow(pv)
  }
  else{
    pv <- apply(pgrid, 2, unique)
    p <- if(is.list(pv)) length(pv[[2]]) else nrow(pv)
  }
  
  omega <- sapply(1:ncol(pgrid), function(i) unique(pgrid[,i])*2 * pi/sidelengths[i])
  omegagrid <- 2*pi*t(t(pgrid)/sidelengths)
  #
  #
  loc <- as.matrix( x$x )
  # "Standardize" coordinates
  if(std){
    loc <- n * t( t( loc/sidelengths ) )
  }
  #
  # Compute:
  CONST <- (1i) * 2 * pi/n
  g <- function(pvec) sum( exp( CONST * loc%*%pvec ) )
  J <- apply(pgrid, 1, g)
  # The constant is chosen mu A = n
  est <- (Re(J)^2 + Im(J)^2) * 2/n
  #
  # Drop (0,0,..) value
  zero <- which( apply(pgrid==0, 1, all) )
  zerov <- est[zero]
  est[ zero ] <- NA
  
  list(v=est, Z=matrix(est, ncol=p ), 
       omega=omega, omegagrid=omegagrid,  
       pv=pv, pgrid=pgrid, 
       std=std, bbox=bbox, zerov=zerov)
}





#' Periodogram as in Mugglestone and Renshaw 2001
#' 
#' @param x point pattern, list of $x and $bbox
#' @param ... ignored
#' 
#' @export

periodogram_MR <- function(x, ...) {
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  bbox <- x$bbox
  dim <- ncol(bbox)
  sidelengths <- bbox_sideLengths(bbox) #apply(bbox, 2, diff)
  n <- nrow(x$x)
  pm <- qm <- round(sqrt(n)/2)
  #
  # The grid of integral values
  pv <- 0:pm
  qv <- setdiff(-qm:qm,0)
  grid <- as.matrix( expand.grid(pv, qv) ) * 2 * pi / n
  #
  #
  loc <- x$x 
  # "Standardize" coordinates
  loc <- sapply(1:dim, function(i) (loc[,i]-bbox[1,i])/sidelengths[i]  )
  #
  # Compute:
  E <- loc%*%t(grid)
  f <- colSums(cos(n*E)^2) + colSums(sin(n*E)^2)
  list(f=f, grid=grid, qv=qv, pv=pv)
}


#' Compute the Discrete Fourier Transform of a point pattern
#'
#' 2D only
#' @param x point pattern
#' @param m maximum steps 
#' @details 
#' The pattern will be rescaled to unit square.
#'   
#' @export
DFT <- function(x, m) {
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  bbox <- x$bbox
  dim <- ncol(bbox)
  sidelengths <- bbox_sideLengths(bbox) #apply(bbox, 2, diff)
  n <- nrow(x$x)
  pm <- qm <- if(missing(m)) round(sqrt(n)/2) else m
  #
  # The grid of integral values
  pv <- 0:pm
  qv <- setdiff(-qm:qm,0)
  grid <- as.matrix( expand.grid(pv, qv) ) * 2 * pi 
  #
  loc <- x$x 
  # "Standardize" coordinates
  loc <- sapply(1:dim, function(i) (loc[,i]-bbox[1,i])/sidelengths[i]  )
  #
  E <- loc%*%t(grid)
  Fo <- colSums(exp(-1i*E))
  out <- list(A=Re(Fo), B=Im(Fo), grid=grid, pv=pv, qv=qv)
  out$f <- out$A^2+out$B^2
  out
}
#' 