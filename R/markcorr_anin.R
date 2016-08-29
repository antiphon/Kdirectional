#' Inhomogeneous anisotropic mark-correlation function, sector version
#' 
#' Estimate a sector-mark correlation function for second order reweighted "inhomogeneous" pattern.
#'
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param marks if x is not marked (x$marks is empty), use these marks
#' @param u unit vector(s) of direction, as row vectors. Default: x and y axis.
#' @param epsilon Central half angle for the directed sector/cone (total angle of the rotation cone is 2*epsilon)
#' @param r radius vector at which to evaluate 
#' @param lambda optional vector of intensity estimates at points
#' @param lambda_h if lambda missing, use this bandwidth in a kernel estimate of lambda(x)
#' @param f test function of the form function(m1, m2) ..., returning a vector of length(m1). default: m1*m2
#' @param r_h smoothing for range dimension, epanechnikov kernel
#' @param stoyan If r_h not given, use r_h=stoyan/lambda^(1/dim). Same as 'stoyan' in spatstat's pcf.
#' @param renormalise See details. 
#' @param bootsize bootstrap size for estimating the normaliser if not given.
#' @param normaliser normalising constant under independent marking. If NULL, estimated with bootstrap.
#' @param border Use translation correction? Default=1, yes. Only for cuboidal windows.
#' @param divisor either "d" or "r". Divide by dist(i,j) ("d") instead of r ("r")?
#' @param ... passed on to e.g. \link{intensity_at_points}
#' @details 
#' 
#' Computes a second order reweighted version of the mark correlation function.
#' 
#' lambda(x) at points can be given, 
#' or else it will be estimated using Epanechnikov kernel smoothing. See 
#' 
#' If 'renormalise=TRUE', we normalise the lambda estimate so that sum(1/lambda(x))=|W|. This corresponds in \code{spatstat}'s \code{Kinhom} to setting 'normpower=2'.
#' 
#' @return 
#' Returns a dataframe.
#' 
#' @useDynLib Kdirectional
#' @export

markcorr_anin <- function(x, marks, u, epsilon, r, lambda=NULL, lambda_h, f = function(a,b) a*b,
                          r_h, stoyan=0.15, renormalise=TRUE, bootsize = 1e5, normaliser = NULL,
                          border=1, divisor = "d", ...) {
  warning("markcorr_anin is experimental.")
  marks <- parse_marks(x, marks)
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
  # 
  # new: handle both divisors in the same function
  fun <- markcor_anin_c
  if(divisor=="r"){
    divisor_i <- 0
  }
  else if(divisor=="d"){
    divisor_i <- 1
  }
  else stop("divisor should 'r' or 'd'")
  # Run
  coord <- x$x
  #
  out <- fun(coord, marks, lambda, bbox, r, u, r_h, epsilon, border, f, divisor_i)
  #
  # The marks under independent marking. Bootstrap average if not pre-known
  mnorm <- if(is.null(normaliser))
    mean(f(sample(marks, bootsize, replace=T), sample(marks, bootsize, replace=T)))
  else normaliser
  #
  # sector pcf estimator normaliser
  norm <- if(dim==2) (4*epsilon) else (4 * pi * (1-cos(epsilon)))
  #
  # in case translation weights are not applied
  if(border==0) norm <- norm * V
  #
  #
  estg <- est <- NULL
  for(i in 1:ncol(out[[1]])) { # per direction
    est <- cbind(est,     (out[[1]][,i]/out[[2]][,i]) / mnorm   ) # all other normalisations factors cancel
    estg <- cbind(estg,   2 * S * out[[2]][,i] / norm  ) # pcf
  }
  #
  # theoretical
  theo <- rep(1, length(r)) 
  #
  # compile a nice output object
  #
  # direction names
  dir_names <- apply(u, 1, function(ui) paste0("(", paste0(ui, collapse=","), ")" ))
  #
  # the mark correlation
  mest <- data.frame(r=r, theo=theo, est)
  names(mest)[] <- c("r", "theo", dir_names)
  rownames(mest) <- NULL
  attr(mest, "epsilon") <- epsilon
  attr(mest, "r_h") <- r_h
  attr(mest, "f") <- f
  class(mest) <- c("markcorr_anin", is(mest))
  #
  # the pair correlation
  gest <- data.frame(r=r, theo=theo, estg)
  names(gest)[] <- c("r", "theo", dir_names)
  rownames(gest) <- NULL
  attr(gest, "epsilon") <- epsilon
  attr(gest, "r_h") <- r_h
  class(gest) <- c("pcf_anin", is(gest))
  # put together
  #
  res <- list(markcorr_anin = mest, pcf_anin = gest)
  class(res) <- c("markcorr_anin", is(res))
  #done
  res
}


#' Plot markcorr_anin object
#' 
#' @param x Output from markcor_anin
#' @param r_scale Plot with x-axis r*r_scale
#' @param rmax plot upto this range
#' @param ylim optional range for y-axis
#' @param legpos legend position
#' @param ... passed on to plot
#' @export

plot.markcorr_anin <- function(x, r_scale=1, rmax, ylim, legpos="topright", ...) {
  if(is.list(x)) x <- x$markcorr_anin
  # cut r
  if(!missing(rmax)) x <- x[x$r<rmax,]
  if(missing(ylim)) ylim <- c(0,2)
  #
  plot(x$r*r_scale, x$theo, col=1, xlab="r", 
       ylab="markcorr_anin", type="l", lty=3, ylim=ylim, ...)
  n <- ncol(x)
  for(i in 3:n){
    lines(x$r*r_scale, x[,i], col=i-1)
  }
  legend(legpos, names(x)[-1], lty=c(3,rep(1,n-2)), col=c(1:(n-1)))
}



