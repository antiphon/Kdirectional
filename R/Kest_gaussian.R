#' Gaussian Weighted Sum of Pairwise Vectors
#' 
#' Estimate K function -type summary for second order reweighted ("inhomogeneous") pattern using a Gaussian kernel sitting in the origin of the Fry-plot.
#'
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param u unit vector(s) of direction, as row vectors. Default: x and y axes, viz. c(1,0) and c(0,1). This gives the major axis of the ellipsoid (direction of largest variance in the Gaussian density).
#' @param kappa The ratio of the secondary axis to the main axis going along 'u'.
#' @param r radius vector at which to evaluate. Corresponds to the radii of the 95\% quantile in x-axis, before rotation in directions u.
#' @param lambda optional vector of intensity estimates at points
#' @param lambda_h if lambda missing, use this bandwidth in a kernel estimate of lambda(x)
#' @param renormalise See details. 
#' @param border Use border correction? Default=1, yes. At the moment no other version available!
#' @param ... passed on to e.g. \link{intensity_at_points}
#' @details 
#'  TODO Gaussian kernel sitting at the fry-space, sum over fry-points.
#' 
#' @return 
#' Returns a dataframe.
#' 
#' @useDynLib Kdirectional
#' @export
Kest_gaussian <- function(x, u, kappa, r, lambda=NULL, 
                          lambda_h, 
                          renormalise=TRUE,  border=1, ...) {
  x <- check_pp(x)
  bbox <- x$bbox
  trans <- !is.bbquad(bbox)
  dim <- bbox_dim(bbox)
  V <- bbox_volume(bbox)
  # directions
  if(missing(u)){ 
    u <- diag(1, dim)
  }
  #
  # make sure unit vectors
  u <- rbind(u)
  u <- t(apply(u, 1, function(ui) ui/c(sqrt(t(ui)%*%ui))))
  #
  # The covariance matrix
  if(missing(kappa)){
    kappa <- 1
  }
  if(kappa > 1 | kappa < 0) stop("kappa should be in range (0, 1]")
  #
  # ranges
  if(missing(r)) {
    sidelengths <- bbox_sideLengths(bbox)
    bl <- min(sidelengths)*0.3
    r <- seq(0, bl, length=50)
  }
  # check intensity 
  if(!missing(lambda)){
    err <- paste("lambda should be a single positive number or a vector of length", nrow(x$x))
    if(!is.vector(lambda)) stop(err)
    
    if(length(lambda) != nrow(x$x)){
      if(length(lambda)!= 1) stop(err)
      lambda <- rep(lambda, nrow(x$x))
    }
  }
  else{
    # estimate lambda
    if(missing(lambda_h)) stop("Need lambda_h to estimate the intensity function")
    lambda <- intensity_at_points(x, bw=lambda_h, ...)
  }
  # Check the lambda's positive
  if(!all(lambda>0)) stop("Check your parameters. Lambda's need to be positive.")
  # if renormalisation of the intensity is in order
  if(renormalise) {
    S <- V/sum(1/lambda)
    normpower <- 2 
    S <- S^normpower
  } else {
    S <- 1
  }
  #
  # 
  # we got everything, let's compute.
  coord <- x$x
  out <- Kest_gaussian_c(coord, lambda, bbox, r, u, kappa, 1)
  out <- 2 * out # double sum
  #
  # scaling if intensity normalised
  out <- S * out
  #
  # compile output
  # direction names
  dir_names <- apply(u, 1, function(ui) paste0("(", paste0(round(ui, 3), collapse=","), ")" ))
  # theoretical under poisson should be just a constant
  theo <- rep(1, length(r))
  #
  Kest <- data.frame(r=r, theo=theo, out)
  names(Kest)[] <- c("r", "theo", dir_names)
  rownames(Kest) <- NULL
  attr(Kest, "kappa") <- kappa
  attr(Kest, "fun_name") <- "Kest_gaussian"
  attr(Kest, "theo_name") <- "CSR"
  #done
  class(Kest) <- c("K_anin", is(Kest))
  Kest
}
