#' Inhomogeneous anisotropic pcf function, cylinder version
#' 
#' Estimate a cylinder-pcf function for second order reweighted ("inhomogeneous") pattern.
#'
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param u unit vector(s) of direction, as row vectors. Default: x and y axis.
#' @param epsilon The cylinder half-width
#' @param r radius vector at which to evaluate the function
#' @param lambda optional vector of intensity estimates at points
#' @param lambda_h if lambda missing, use this bandwidth in a kernel estimate of lambda(x)
#' @param r_h smoothing for range dimension, epanechnikov kernel
#' @param stoyan If r_h not given, use r_h=stoyan/lambda^(1/dim). Same as 'stoyan' in spatstat's pcf.
#' @param renormalise See details. 
#' @param border Use translation correction? Default=1, yes. Only for cuboidal windows.
#' @param ... passed on to e.g. \link{intensity_at_points}
#' @details 
#' 
#' Computes a second order reweighted version of the cylinder-pcf, defined as the function to integrate in range over [0,R] to get the cylinder-K(R) function.
#' 
#' Lambda(x) at points can be given, 
#' or else it will be estimated using Epanechnikov kernel smoothing. See 
#' 
#' If 'renormalise=TRUE', we normalise the lambda estimate so that sum(1/lambda(x))=|W|. This corresponds in \code{spatstat}'s \code{Kinhom} to setting 'normpower=2'.
#' 
#' @return 
#' Returns a dataframe.
#' 
#' @useDynLib Kdirectional
#' @export

pcf_anin_cylinder <- function(x, u, epsilon, r, lambda=NULL, lambda_h, r_h, 
                              stoyan=0.15,
                     renormalise=TRUE,  border=1, ...) {
  x <- check_pp(x)
  bbox <- x$bbox
  if(is.bbquad(bbox)) stop("bbquad window not yet supported.")
  dim <- ncol(bbox)
  V <- bbox_volume(bbox)
  # directions
  if(missing(u)){ 
    u <- diag(c(1), dim)
  }
  #
  # make sure unit vectors
  u <- rbind(u)
  u <- t(apply(u, 1, function(ui) ui/sqrt(t(ui)%*%ui)))
  #
  # central half-angle
  if(missing(epsilon)){
    epsilon <- 0.1 
  }
  if(epsilon<= 0) stop("Cylinder radius epsilon should be >0")
  #
  # ranges
  if(missing(r)) {
    sidelengths <- apply(bbox, 2, diff)
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
  
  
  if(missing(r_h)) {
    lambda0 <- nrow(x$x)/V # mean lambda
    r_h <- stoyan/lambda0^(1/dim)
  }
  
  # if renormalisation of the intensity is in order
  if(renormalise) {
    S <- V/sum(1/lambda)
    normpower <- 2 
    S <- S^normpower
  } else {
    S<-1
  }
  # new, do it in one function
  # check divisor
  fun <- pcf_anin_cylindrical_c
  div <- 1
  # Run
  coord <- x$x
  out <- fun(coord, lambda, bbox, r, r_h, u, epsilon, border)
  #
  # scaling
  #
  out <- S * out # double sum
  # in case translation weights are not applied
  if(border==0) out <- out/V
  #
  # scale
  norm <- pi^((dim-1)/2)/gamma((dim-1)/2+1) * epsilon^(dim-1) * div
  out <- out/norm
  #
  # compile output
  # direction names
  dir_names <- apply(u, 1, function(ui) paste0("(", paste0(ui, collapse=","), ")" ))
  # theoretical
  theo <- rep(1, length(r)) 
  #
  gest <- data.frame(r=r, theo=theo, out)
  names(gest)[] <- c("r", "theo", dir_names)
  rownames(gest) <- NULL
  attr(gest, "epsilon") <- epsilon
  attr(gest, "r_h") <- r_h
  attr(gest, "fname") <- "pcf_anin_cylinder"
  #done
  class(gest) <- c("pcf_anin", is(gest))
  gest
}

