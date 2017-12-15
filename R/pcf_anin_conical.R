#' Inhomogeneous anisotropic pcf function, sector version
#' 
#' Estimate a sector/cone pcf function for second order reweighted ("inhomogeneous") pattern.
#'
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param u unit vector(s) of direction, as row vectors. Default: x and y axis.
#' @param epsilon Central half angle for the directed sector/cone (total angle of the rotation cone is 2*epsilon)
#' @param r radius vector at which to evaluate the function
#' @param lambda optional vector of intensity estimates at points
#' @param lambda_h if lambda missing, use this bandwidth in a kernel estimate of lambda(x)
#' @param r_h smoothing for range dimension, epanechnikov kernel
#' @param stoyan If r_h not given, use r_h=stoyan/lambda^(1/dim). Same as 'stoyan' in spatstat's pcf.
#' @param renormalise See details. 
#' @param border Use translation correction? Default=1, yes. Only for cuboidal windows.
#' @param divisor See spatstat's pcf.ppp for this.
#' @param ... passed on to e.g. \link{intensity_at_points}
#' @details 
#' 
#' Computes a second order reweighted version of the sector-pcf. The sector pcf differs from the true anisotropic pcf
#' by assuming that the anisotropic pcf is constant over the small arc/cap of the sector, thus averaging over that
#' data-area and providing more stable estimates. 
#' 
#' Lambda(x) at points can be given, 
#' or else it will be estimated using kernel smoothing. See \link{intensity_at_points}.
#' 
#' If 'renormalise=TRUE', we normalise the lambda estimate so that sum(1/lambda(x))=|W|. This corresponds in \code{spatstat}'s \code{Kinhom} to setting 'normpower=2'.
#' 
#' @return 
#' Returns a dataframe.
#' 
#' @seealso \code{\link{pcf_anin_cylinder}}
#' 
#' @useDynLib Kdirectional
#' @export

pcf_anin <- pcf_anin_conical <- function(x, u, epsilon, r, lambda=NULL, lambda_h, r_h, stoyan=0.15,
                      renormalise=TRUE,  border=1, divisor = "d", ...) {
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
  u <- t(apply(u, 1, function(ui) ui/c(sqrt(t(ui)%*%ui))))
  #
  # central half-angle
  if(missing(epsilon)){
    epsilon <- pi/4 
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
  if(divisor=="r"){
    divisor_i <- 0
    div <- 1 
  }
  else if(divisor=="d"){
    divisor_i <- 1  
    div <- 1 
  }
  else stop("divisor should 'r' or 'd'")
  # Run
  coord <- x$x
  
  fun <- pcf_anin_conical_c # handles both divisors inside
  
  out <- fun(coord, lambda, bbox, r, r_h, u, epsilon, border, divisor_i)
  #
  # scaling
  #
  out <- 2  * S * out # double sum
  # in case translation weights are not applied
  if(border==0) out <- out/V
  #
  # scale
  norm <- if(dim==2) (4*epsilon*div) else (4 * div * pi * (1-cos(epsilon)))
  if(dim > 3) warning("Dimension beyond 3, normalising constant unknown.")
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
  attr(gest, "fname") <- "pcf_anin_conical"
  #done
  class(gest) <- c("pcf_anin", is(gest))
  gest
}


#' Plot pcf_anin object
#' 
#' @param x Output from pcf_anin or pcf_anin_cylinder
#' @param r_scale Plot with x-axis r*r_scale
#' @param rmax plot upto this range
#' @param ylim optional range for y-axis
#' @param legpos legend position
#' @param ... passed on to plot
#' @export

plot.pcf_anin <- function(x, r_scale=1, rmax, ylim, legpos="topright", lwd = lwd, ...) {
  # cut r
  if(!missing(rmax)) x <- x[x$r<rmax,]
  if(missing(ylim)) ylim <- c(0,2)
  #
  fname <- attr(x, "fname")
  plot(x$r*r_scale, x$theo, col=1, xlab="r", 
       ylab=if(is.null(fname)) "pcf_anin" else fname, type="l", lty=3, ylim=ylim, lwd = lwd, ...)
  n <- ncol(x)
  for(i in 3:n){
    lines(x$r*r_scale, x[,i], col=i-1, lwd = lwd)
  }
  legend(legpos, names(x)[-1], lty=c(3,rep(1,n-2)), col=c(1:(n-1)), lwd = lwd)
}






