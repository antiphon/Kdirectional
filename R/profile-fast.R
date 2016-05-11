#' Anisotropy profile for given diagonal transformations
#' 
#' Compares the directed sector-K functions along main axes.
#' 
#' @param x point pattern
#' @param grid p x dim matrix giving the diagonals of the candidate transformations
#' @param r vector of range values to sum over ('integration nodes')
#' @param verb verbose?
#' @param power Power of the differences, absolute value for power=1.
#' @param eps Sector angle
#' @param border Border correction? Should use this unless the data is large.
#' @param antipodal Use antipodally symmetric estimation of sector-K?
#' 
#' @details
#' For each v in grid, transform x -> T(v)x and summarize anisotropy by
#' absolute integral difference in x-y-z-directed K-functions.
#' 
#' The transformations T(v) are given by diagonal vector v, rows of parameter \code{grid}.
#' 
#' We assume that no rotation is present, or that rotation is already corrected for, so that
#' squared difference sums of the directed sector-K functions along main axes should indicate
#' the optimal v.
#' 
#' @export
#' 

anisotropy_profile_fast <- function(x, grid, r, verb=FALSE, eps=pi/4, border=TRUE, antipodal=TRUE, power=2, ...) {
  x <- check_pp(x)
  
  if(!is.bbquad(x$bbox)) x$bbox <- bbox2bbquad(x$bbox)
  
  dim <- ncol(x$x)
  if(is.null(ncol(grid))) stop("grid should be matrix.")
  if(dim != ncol(grid)) stop(paste("grid should a p x", dim, " matrix."))
  
  if(missing(r)) r <- seq(0, min(bbox_sideLengths(x$bbox))*0.3, length=50)
  
  cat2 <- if(verb) cat else function(...)NULL

  # if no active border correction
  if(!border){
    fry <- fry_points(x, border = TRUE)
    Fxyz <- fry$fry_r * fry$fry_units
  }
  
  # Inverse transform the fry points, need to recompute 
  # since border distances change.
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
  scaling <- 1/(nrow(x$x)/vol) # for scaling the K function
  
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
      val <- anisotropy_abs_fast(Ff, r, eps, antipodal, power)
      vals <- append(vals, list(val))
      # gather
      profile <- c(profile, val$stat * scaling)
    }
    cat2(k<-k+1, "/", length(gridl), "   \r")
  }
  cat2("\nDone.")
  #profile <- data.frame(grid=grid, anisotropy=profile)
  out <- list(grid=grid, profile=profile, 
              raw=vals, r=r, power=power,
              eps=eps, border=border, antipodal=antipodal)
  class(out) <- "anisotropyprofile"
  out
}

#' print for anisotropy profile
#' 
#' @export
print.anisotropyprofile <- function(x, ...) {
  print( summary(x,...) )
}

#' Summary for anisotropy profile
#' 
#' @export
summary.anisotropyprofile <- function(x, ...) {
  y <- x
  tab <- x$profile
  i <- which.min(tab)
  x$estimate <- diag(x$grid[i,])
  class(x) <- "anisotropyprofile.summary"
  x
}

#' Print method for "anisotropyprofile" summary
#' 
#' @export
print.anisotropyprofile.summary <- function(x, ...){
  cat("Anisotropy profile\n")
  cat("Range of integration: [", min(x$r), ",", max(x$r), "]\n")
  cat("Difference power:", x$power, "\n")
  cat("Optimal transformation: diag(", diag(x$estimate), ")" )
}

#' Plot method for "anisotropyprofile"
#' 
#' @export
plot.anisotropyprofile <- function(x, ..., 
                                   xlab = "first grid component",
                                   ylab = "anisotropy profile",
                                   scale = TRUE) {
  d <- ncol(x$grid)
  if(scale) x$profile <- (x$profile-min(x$profile))/diff(range(x$profile))
  if(d == 2) {
    plot(x$grid[,1], x$profile, xlab=xlab, ylab=ylab , ...)
  }
  else{ # assume 3D
    stop("Plot not implemented for 3D")
  }
}  

#' Lines method for "anisotropyprofile"
#' 
#' @export
lines.anisotropyprofile <- function(x, ..., scale=TRUE) {
  d <- ncol(x$grid)
  if(scale) x$profile <- (x$profile-min(x$profile))/diff(range(x$profile))
  if(d == 2) {
    lines(x$grid[,1], x$profile, ...)
  }
  else{ # assume 3D
    stop("Lines not implemented for 3D")
  }
}  


#####################################################################
#' Fast absolute anisotropy  statistic
#'
#' @param f Fry points, border corrected
#' @param r range of integration
#' @param eps width of cone angle
#' @param antipodal If TRUE, antipodes are equated
#'
anisotropy_abs_fast <- function(f, r, eps, antipodal, power){
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
  metr <- if(power==1) abs else function(x) x^power
  if(is.null(nrow(inside))){
    warning("No data to bin.")
    counts <- matrix(rep(0, dim*length(r)), ncol=dim)
  }
  else {
    counts <- apply(inside, 2, function(ins) sapply(r, function(rr) sum(Fr[ins]<rr) ))
  }
  
  if(dim==2){
    stat <- sum(metr(counts[,1]-counts[,2]))
  }
  else{
    stat <- sum( metr(counts[,1]-counts[,2]) + metr(counts[,1]-counts[,3]) + metr(counts[,2]-counts[,3])     )
  }
  # add dr
  stat <- stat * mean(diff(r))
  out <- list(statistic=stat, raw=counts)
  out
}
