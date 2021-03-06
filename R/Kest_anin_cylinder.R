#' Anisotropic and Inhomogeneous K function, cylinder version
#' 
#' Estimate a cylinder-K function for second order reweighted ("inhomogeneous") pattern.
#'
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param u unit vector(s) of direction, as row vectors. Default: x and y axes, viz. c(1,0) and c(0,1).
#' @param epsilon The cylinder half-width. Will be extended to the length of r, so can be given per r.
#' @param r radius vector at which to evaluate K
#' @param lambda optional vector of intensity estimates at points
#' @param lambda_h if lambda missing, use this bandwidth in a kernel estimate of lambda(x)
#' @param renormalise See details. 
#' @param border Use border correction? Default=1, yes. 
#' @param aspect Instead of using a fixed halfwidth (epsilon) take the halfwidth to be 'aspect * r/2' (so an increasing vector). Default : 1/3
#' @param ... passed on to e.g. \link{intensity_at_points}
#' @details 
#' 
#' Computes a second order reweighted version of the cylinder-K. In short, we count how many pairs of points in the pattern 
#' has both a) their difference vector inside a cylinder with major-axial direction 'u' and radius epsilon,  and 
#' b) difference vector length less than range r. Usually r is a vector and the output is then a vector as well.
#' 
#' Note the default behaviour is to use a fixed aspect ratio cylinder with aspect = 2*epsilon/range = 1/3.
#' 
#' An estimate of the intensity Lambda(x) at points can be given ('lambda'). If it is a single value, the pattern is assumed to be homogeneous. 
#' If it is a vector the same length as there are points, the pattern is taken to be second-order stationary. In this case the 
#' the sum over the pairs (i,j) is weighted with 1/(lambda[i]*lambda[j]). If 'lambda' is missing, 'lambda_h', a single positive number, 
#' should be given, which is then used for estimating the non-constant Lambda(x) via Epanechnikov kernel smoothing (see \link{intensity_at_points}).
#' If 'renormalise=TRUE', we normalise the intensity estimate so that sum(1/lambda(x))=|W|. This corresponds in \code{spatstat}'s \code{Kinhom} to setting 'normpower=2'.
#' 
#' About border correction: If x$bbox is a a simple bounding box, the algorithm uses the translation corrected weighting 1/area(Wx intersect Wy) with Wx=W+x. If x$bbox is a bbquad-object, for example rotated polygon, the algorithm uses simple minus border correction.
#' 
#' 
#' @return 
#' Returns a dataframe.
#' 
#' @useDynLib Kdirectional
#' @export

Kest_anin_cylinder <- function(x, u, epsilon, r, lambda=NULL, lambda_h, 
                      renormalise=TRUE,  border=1, aspect = 1/3, ...) {
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
  u <- t(apply(u, 1, function(ui) ui/c(sqrt(t(ui)%*%ui) )))
  #
  #
  # ranges
  if(missing(r)) {
    sidelengths <- bbox_sideLengths(bbox)
    bl <- min(sidelengths)*0.3
    r <- seq(0, bl, length=50)
  }
  
  # central half-angle
  if(!is.null(aspect)) {
    epsilon <- aspect * r / 2
  }
  else{
    if(missing(epsilon)){
      stop("cylinder half-width epsilon missing, and no default available.")
    }
    if(length(epsilon) == 1) epsilon <- rep(epsilon, length(r)) 
    else if( length(epsilon) == length(r)) {}
    else stop("epsilon should be length 1 or length(r).") 
  }
  if(any(epsilon < 0)) stop("epsilon should >0")
  #
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
    S<-1
  }
  
  # 
  # we got everything, let's compute.
  coord <- x$x
  if(trans){
    out <- Kest_anin_cylinder_c(coord, lambda, bbox, r, u, epsilon, border)
    S <- 2 * S # double sum
  }
  if(!trans){ # minus border correction
    bd <- bbox_distance(coord, bbox)
    bbox0 <- bbquad2bbox(bbox)
    out <- Kest_anin_cylinder_border_c(coord, lambda, bbox0, bd, r, u, epsilon, border)
    # need to correct due to minus sampling
    if(border){
      w <- sapply(r, function(r) sum(bd > r))
      out <- out*nrow(coord)/w
    }
  }
  # scaling
  #
  out <- S * out
  # in case translation weights are not applied
  if(border==0 |  !trans ) out <- out/V 
  #
  # compile output
  # direction names
  dir_names <- apply(u, 1, function(ui) paste0("(", paste0(round(ui, 3), collapse=","), ")" ))
  # theoretical
  theo <- if(dim==2) (4 * epsilon * r) else (2 * epsilon^2 * pi * r)
  #
  Kest <- data.frame(r=r, theo=theo, out)
  names(Kest)[] <- c("r", "theo", dir_names)
  rownames(Kest) <- NULL
  attr(Kest, "epsilon") <- epsilon
  attr(Kest, "aspect")  <- aspect
  attr(Kest, "fun_name") <- "Kest_anin_cylinder"
  attr(Kest, "theo_name") <- "CSR"
  
  #done
  class(Kest) <- c("K_anin", is(Kest))
  Kest
}
