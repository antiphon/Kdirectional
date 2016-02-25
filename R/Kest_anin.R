#' Inhomogeneous anisotropic K function
#' 
#' Estimate a Sector-K function for second order reweighted ("inhomogeneous") pattern.
#'
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param u unit vector(s) of direction, as row vectors. Default: x and y axis.
#' @param epsilon Central half angle for the directed sector/cone (total angle of the rotation cone is 2*epsilon)
#' @param r radius vector at which to evaluate K
#' @param lambda optional vector of intensity estimates at points
#' @param lambda_h if lambda missing, use this bandwidth in a kernel estimate of lambda(x)
#' @param renormalise See details. 
#' @param border Use translation correction? Default=1, yes. Only for cuboidal windows.
#' @param ... passed on to e.g. \link{intensity_at_points}
#' @details 
#' 
#' Computes a second order reweighted version of the Sector-K.
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

Kest_anin <- function(x, u, epsilon, r, lambda=NULL, lambda_h, 
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
   epsilon <- pi/2 
  }
  if(abs(epsilon)>pi/2) stop("epsilon should be in range [0, pi/2]")
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
  out <- Kest_anin_c(coord, lambda, bbox, r, u, epsilon, border)
  # scaling
  #
  out <- 2  * S * out # double sum
  # in case translation weights are not applied
  if(border==0) out <- out/V 
  #
  # compile output
  # direction names
  dir_names <- apply(u, 1, function(ui) paste0("(", paste0(ui, collapse=","), ")" ))
  # theoretical
  theo <- if(dim==2) (2*epsilon*r^2) else (4/3 * r^3 * pi * (1-cos(epsilon)))
  #
  Kest <- data.frame(r=r, theo=theo, out)
  names(Kest)[] <- c("r", "theo", dir_names)
  rownames(Kest) <- NULL
  attr(Kest, "epsilon") <- epsilon
  #done
  class(Kest) <- c("K_anin", is(Kest))
  Kest
}


#' Plot Kest_anin object
#' 
#' @param x Output from Kest_anin
#' @param r_scale Plot with x-axis r*r_scale
#' @param rmax plot upto this range
#' @param legpos legend position
#' @param ... passed on to plot
#' @export

plot.K_anin <- function(x, r_scale=1, rmax, legpos="topleft", ...) {
  # cut r
  if(!missing(rmax)) x <- x[x$r<rmax,]
  #
  plot(x$r*r_scale, x$theo, col=1, xlab="r", 
       ylab="Kest_anin", type="l", lty=3, ...)
  n <- ncol(x)
  for(i in 3:n){
    lines(x$r*r_scale, x[,i], col=i-1)
  }
  legend(legpos, names(x)[-1], lty=c(3,rep(1,n-2)), col=c(1:(n-1)))
}


