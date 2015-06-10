#' Anisotropy profile for given diagonal transformations
#' 
#' @param x point pattern, list(x=coords, bbox=bounding-box)
#' @param grid p x dim matrix giving the diagonals of the candidate transformations
#' 
#' @details
#' For each v in grid, transform x -> T(v)x and summarize anisotropy by
#' absolute integral difference in x-y-z-directed K-functions.
#' 
#' The transformations T(v) are given by diagonal vector v, rows of parameter \code{grid}.
#' 
#' We assume that no rotation is present, or that rotation is already corrected for. 
#' 
#' @export
#' 

anisotropy_profile_fast <- function(x, grid, r, verb=FALSE, eps=pi/4, border=TRUE, antipodal=TRUE, ...) {
  x <- check_pp(x)
  dim <- ncol(x$x)
  if(is.null(ncol(grid))) stop("grid should be matrix.")
  if(dim != ncol(grid)) stop(paste("grid should a p x", dim, " matrix."))
  
  if(missing(r)) r <- seq(0, min(bbox_sideLengths(x$bbox))*0.3, length=50)
  
  cat2 <- if(verb) cat else function(...)NULL

  fry <- fry_points(x, border = border)
  Fxyz <- fry$fry_r * fry$fry_units
  
  
  #' inverse transform the fry points
  Theta <- function(ce) {
    A <- diag(1/ce)
    x0 <- x$x%*%t(A)
    if(border){
      z <- list(x=x0, bbox=bbox_affine(x$bbox, A))
      fry <- fry_points(z, border = TRUE)
      fry$fry_r * fry$fry_units
    } else Fxyz%*%t(A)
  }
  
  k<-0
  vals <- list()
  gridl <- split(grid, 1:nrow(grid))
  profile <- NULL
  vals <- list()
  
  vol <- bbox_volume(x$bbox) 
  scaling <- 1/(nrow(x$x)/vol)
  
  for(ce in gridl) {
    # deform fry points
    Ff <- Theta(ce)
    N <- nrow(Ff)
    if(is.null(N)){
      warning(paste0("No valid fry points left after scaling:", paste0(ce, collapse=",")))
      profile <- c(profile, NA)
    }
    else{
      # compute anisotropy using fry points:
      val <- anisotropy_abs_fast(Ff, r, eps, antipodal)
      vals <- append(vals, list(val))
      # gather
      profile <- c(profile, val$stat * scaling)
    }
    cat2(k<-k+1, "/", length(gridl), "   \r")
  }
  cat2("\nDone.")
  profile <- data.frame(grid=grid, anisotropy=profile)
  
  list(profile=profile, raw=vals)
}



#####################################################################
#' Fast absolute anisotropy  statistic
#'
#' @param f Fry points, border corrected
#' @param r range of integration
#' @param eps width of cone angle
#' @param antipodal If TRUE, antipodes are equated
#'
#' @export
anisotropy_abs_fast <- function(f, r, eps, antipodal){
  Fr <- sqrt(rowSums(f^2))
  # drop too longs
  o <- Fr < max(r)
  f <- f[o,]
  Fr <- Fr[o]
  #
  Fu <- f/Fr
  dim <- ncol(f)
  
  # compute unnormalized K's
  ang <- acos(Fu %*% diag(1,dim,dim))
  if(antipodal) ang <- pmin(ang, pi-ang)
  inside <- ang < eps
  if(is.null(nrow(inside))){
    warning("No data to bin.")
    counts <- matrix(rep(0, dim*length(r)), ncol=dim)
  }
  else {
    counts <- apply(inside, 2, function(ins) sapply(r, function(rr) sum(Fr[ins]<rr) ))
  }
  
  if(dim==2){
    stat <- sum(abs(counts[,1]-counts[,2]))
  }
  else{
    stat <- sum( abs(counts[,1]-counts[,2])+abs(counts[,1]-counts[,3])+abs(counts[,2]-counts[,3])     )
  }
  # add dr
  stat <- stat * mean(diff(r))
  
  list(statistic=stat, raw=counts)
}




