#' Anisotropic and Inhomogeneous T function
#' 
#' Estimate a Sector-T function for second order reweighted ("inhomogeneous") pattern.
#'
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param u unit vector(s) of direction, as row vectors. Default: x and y axes, viz. c(1,0) and c(0,1).
#' @param epsilon Central half angle for the directed sector/cone (total angle of the rotation cone is 2*epsilon). Default: pi/4.
#' @param r radius vector at which to evaluate K
#' @param lambda optional vector of intensity estimates at points
#' @param lambda_h if lambda missing, use this bandwidth in a kernel estimate of lambda(x)
#' @param renormalise See details. 
#' @param border Use border correction? Default=1, yes. 
#' @param ... passed on to e.g. \link{intensity_at_points}
#' @details 
#' 
#' Computes a second order reweighted version of the Sector-T. In short, we count how many triplets of points in the pattern 
#' have both a) their difference vector's angle less than 'epsilon' radians from direction 'u' and 
#' b) difference vector lengths less than range r. Usually r is a vector and the output is then a vector as well.
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

Test_anin <- function(x, u, epsilon, r, lambda=NULL, lambda_h, 
                      renormalise=TRUE,  border=1, ...) {
  stop("Not implemented.")
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
  # central half-angle
  if(missing(epsilon)){
    epsilon <- pi/4
  }
  if(abs(epsilon)>pi/2) stop("epsilon should be in range [0, pi/2]")
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
    S<-1
  }
  
  # 
  # we got everything, let's compute.
  coord <- x$x
  if(trans){
    stop("Not implemented.")
    out <- Kest_anin_c(coord, lambda, bbox, r, u, epsilon, border)
    out <- out * 2 # double sum
  }
  if(!trans){ # minus border correction
    bd <- bbox_distance(coord, bbox)
    bbox0 <- bbquad2bbox(bbox)
    xout <- Kest_anin_border_c(coord, lambda, bbox0, bd, r, u, epsilon, border)
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
  theo <- if(dim==2) (2 * epsilon * r^2) else (4/3 * r^3 * pi * (1-cos(epsilon)))
  #
  Kest <- data.frame(r=r, theo=theo, out)
  names(Kest)[] <- c("r", "theo", dir_names)
  rownames(Kest) <- NULL
  attr(Kest, "epsilon") <- epsilon
  attr(Kest, "fun_name") <- "Kest_anin"
  attr(Kest, "theo_name") <- "CSR"
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

plot.K_anin <- function(x, r_scale=1, rmax, legpos="topleft", lwd = 1, ...) {
  # cut r
  if(!missing(rmax)) x <- x[x$r<rmax,]
  #
  plot(x$r*r_scale, x$theo, col=1, xlab="r", 
       ylab=attr(x, "fun_name"), type="l", lty=3, lwd = lwd, ...)
  n <- ncol(x)
  for(i in 3:n){
    lines(x$r*r_scale, x[,i], col=i-1, lwd = lwd)
  }
  nam <- names(x)[-1]
  tn <- attr(x, "theo_name")
  if(!is.null(tn)) nam[1] <- tn
  legend(legpos, nam, lty=c(3,rep(1,n-2)), col=c(1:(n-1)), lwd = lwd)
}



