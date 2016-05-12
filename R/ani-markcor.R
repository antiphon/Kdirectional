#' Anisotropic mark correlation function, direction vector formulation
#' 
#' Estimate the anisotropic mark correlation function (2d and 3d), as defined in Stoyan 1991, f. 5.4. (without the additional 'a'-marks)
#' 
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param marks if x is not marked (x$marks is empty), use these marks
#' @param r Evaluate at these lengths.
#' @param f test function of the form function(m1, m2) ..., returning a vector of length(m1)
#' @param directions Matrix of directions, in unit vectors, one direction per row. Default: along axes.
#' @param bw Bandwidths of epanechnicov kernels. Vector of two values, one for ranges and one for angles.
#' @param adjust If bw not given, use bw=adjust[1]*0.15/lambda^(1/dim) for range and bw=adjust[2]*0.15*pi for angle
#' @param correction "none" or "translation". Translation only for rectangular box.
#' @param bootsize bootstrap size for estimating the normaliser if not given.
#' @param divisor either "d" or "r". Divide by dist(i,j) ("d") instead of r ("r")?
#' @param normaliser normalising constant under independent marking. If NULL, estimated with bootstrap.
#' 
#' @details
#' 
#' 
#' 
#' Default bandwidth: bw=adjust*0.15/lambda^(1/dim) for range and bw=adjust*0.15*pi for angle.
#' 
#' 
#' 
#' @useDynLib Kdirectional
#' @import rgl sphere
#' @export
markcorr_anisotropic <- function(x, marks=NULL, r, f = function(a,b) a*b, 
                               directions,
                               bw, adjust=c(1,1), 
                               correction="translation", 
                               bootsize=1e5,
                               divisor = "d",
                               normaliser = NULL
                               ) {
  
  marks <- parse_marks(x, marks)
  x <- check_pp(x)
  bbox <- as.matrix(x$bbox)
  dim <- ncol(bbox)
  sidelengths <- bbox_sideLengths(bbox)
  lambda <- nrow(x$x)/bbox_volume(bbox)
  #
  #
  if(missing(r)){ 
    b <- min(sidelengths)*0.2
    r <- seq(0, b, length=50)
  }
  # Check directions
  if(missing(directions)){
    directions <- diag(rep(1, dim))
  }
  else{# directions given. make sure unit vectors.
    directions <- rbind(directions)
    if(ncol(directions)!=dim) stop(paste("'directions' should be in dimension", dim ))
    dl <- sqrt(rowSums(directions^2))
    directions <- directions/dl
  }
  #
  #
  #
  # check smoothing parameters
  if(missing(bw)) {
    bw <- adjust * c(0.15/lambda^(1/dim), 0.15*pi)
  }
  if(length(bw) < 2) stop("'bw' should be of length 2.")
  #
  #
  #
  # Edge correction
  correction_i <- pmatch(correction, c("none", "translation"))
  #
  # start:
  xc <- as.matrix(x$x)
  
  if(divisor=="r"){
    fun <- anisotropic_markcor_c
    div <- r^(dim-1) * 0.5^(dim-2)
  }
  else if(divisor=="d"){
    fun <- anisotropic_markcor_c_d
    div <- 1  * 0.5^(dim-2)
  }
  else stop("divisor should 'r' or 'd'")
  # Run
  res <- fun(xc, marks, bbox, r, directions, bw[1], bw[2], f)
  #
  # pcf scale
  scaleg <- div * lambda^2
  #
  # Normalizing constant, use bootstrap to estimate:
  mnorm <- if(is.null(normaliser))
              mean(f(sample(marks, bootsize, replace=T), sample(marks, bootsize, replace=T)))
           else normaliser
  #
  # The estimates:
  estg <- est <- NULL
  for(i in 1:ncol(res[[1]])) { # per direction
    est <- cbind(est,   (res[[1]][,i]/res[[2]][,i]) / mnorm   )
    estg <- cbind(estg,  res[[2]][,i]/scaleg   )
  }
  #
  # done
  list(r=r, directions=directions, mcor=est, pcf=estg, bw=bw)
}


