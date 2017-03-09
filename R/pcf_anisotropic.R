#' Anisotropic pair correlation function
#' 
#' Estimate the anisotropic pair correlation function (2d and 3d), as defined in f 1991, f. 5.2-5.3.
#' 
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param r radius vector at which to estimate
#' @param u direction to be converted to angles. Use \code{\link{angle_2_unit}} to transform from angles.
#' @param h half-widths of epanechnicov kernels, vector of two values, one for ranges and and one for angles.
#' @param stoyan If h not given, use h=c( stoyan/lambda^(1/dim), stoyan*pi) (cf. \code{pcf.ppp}-function in \code{spatstat}-package).
#' @param lambda optional vector of intensity estimates at points
#' @param lambda_h if lambda missing and lambda_h is given, use this bandwidth in a kernel estimate of lambda(x). Otherwise lambda is set to constant.
#' @param border Use translation correction? Default=1, yes. Only for cuboidal windows.
#' @param renormalise Scale lambda to align with Campbell formula. 
#' 
#' @details 
#' The antipode symmetry makes it necessary to compute only on one half of the circle/sphere.
#' 
#' @useDynLib Kdirectional
#' @export

pcf_anisotropic <- function(x, r, u, 
                            h, stoyan=0.15,
                            lambda, lambda_h, renormalise=TRUE, 
                            border=1, divisor = "d", ...) {
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
  # ranges
  if(missing(r)) {
    sidelengths <- apply(bbox, 2, diff)
    bl <- min(sidelengths)*0.3
    r <- seq(0, bl, length=50)
  }
  # check intensity 
  err <- paste("lambda should be a single positive number or a vector of length", nrow(x$x))
  if(!missing(lambda)){
    if(!is.vector(lambda)) stop(err)
  }
  else{
    # estimate lambda
    if(missing(lambda_h)) lambda <- nrow(x$x)/V # constant
    else lambda <- intensity_at_points(x, bw=lambda_h, ...) # not constant
  }
  if(length(lambda) != nrow(x$x)){
    if(length(lambda)!= 1) stop(err)
    lambda <- rep(lambda, nrow(x$x))
  }
  
  # Check the lambda's positive
  if(!all(lambda>0)) stop("Check your parameters. Lambda's need to be positive.")
  
  
  if(missing(h)) {
    lambda0 <- nrow(x$x)/V # mean lambda
    r_h <- stoyan * c(1/lambda0^(1/dim))
    a_h <- stoyan * pi
  }
  else {
    if(length(h)!=2) stop("h should be vector of length 2.")
    r_h <- h[1]
    a_h <- h[2]
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
  
  out <- fun(coord, lambda, bbox, r, r_h, u, a_h, border, divisor_i, ang_kernel = 1)
  #
  # scaling
  #
  out <- 2  * S * out # double sum
  # in case translation weights are not applied
  if(border==0) out <- out/V
  #
  # scale
  #browser()
  norm <- if(dim==2){
    out <- out/2
  }else if(dim == 3){
    out <- out/(2 * a_h)  #t(t(out)/sin( unit_2_angle(u)[,2] ))/2
    warning("Normalisation biased.")
  }
  else warning("Dimension beyond 3, normalising constant unknown.")
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
  attr(gest, "bw_r") <- r_h
  attr(gest, "bw_ang") <- a_h
  attr(gest, "fname") <- "pcf_anisotropic"
  #done
  class(gest) <- c("pcf_anin", is(gest))
  gest
}



#' Anisotropic pair correlation function, old version
#' 
#' Estimate the anisotropic pair correlation function (2d and 3d), as defined in f 1991, f. 5.2-5.3.
#' 
#' @param x pp, list with $x~coordinates $bbox~bounding box
#' @param r radius vector at which to estimate
#' @param theta vector or list of vectors for for angles, see details.
#' @param h widths of epanechnicov kernels, vector of two values, one for ranges and and one for angles.
#' @param stoyan If h not given, use h=c( stoyan/lambda^(1/dim), stoyan*pi) (cf. \code{pcf.ppp}-function in \code{spatstat}-package).
#' @param correction "none" or translation. Translation only for rectangle box.
#' @param n_dir Angle grid resolution (in case theta not given)
#' 
#' The antipode symmetry makes it necessary to compute only on one half of the circle/sphere.
#' 
#' @useDynLib Kdirectional
#' @export

pcf_anisotropic_old <- function(x, r, theta, h, stoyan=0.15, correction="translation", n_dir=7) {
  x <- check_pp(x)
  bbox <- as.matrix(x$bbox)
  dim <- ncol(bbox)
  sidelengths <- apply(bbox, 2, diff)
  lambda <- nrow(x$x)/bbox_volume(bbox)
  
  if(missing(theta)){
    dirs <- check_directions(n_dir=n_dir ,dim=dim)
    theta <- dirs$theta
    if(dim==3) theta <- append(theta, theta)
  }
  if(!is.list(theta)) stop("Theta should be a list.")
  if(any(sapply(theta, abs)>pi)) stop("Theta should be in range [0, pi]")
  if(missing(r)) {
    b <- min(sidelengths)*0.3
    r <- seq(0, b, length=50)
  }
  if(missing(h)) {
    f <- stoyan
    h <- c( f/lambda^(1/dim), f*pi )  
  }
  if(length(h) < 2) h <- rep(h, 2)
  # correction
  correction_i <- pmatch(correction, c("none", "translation"))
  # start:
  xc <- as.matrix(x$x)
  #
  units <- theta_2_unit(theta)
  #
  directions <- NULL 
  for(ri in r) directions <- rbind(directions, ri*units)
  #
  res <- c_anisotropic_unit_pcf(xc, directions, h, bbox, correction_i)
  #  correction
  grid <- expand.grid(append(list(r=r), theta))
  
  rho <- res[[1]]
  # pcf
  g <- rho/lambda^2
  g_m <- matrix(g, byrow=T, nrow=length(r))
  #
  res <- list(est=g_m, r=r, theta=theta, directions=directions, correction=correction, h=h, rho=res[[1]], dim=dim, counts=res[[2]])
  class(res) <- "pcf_anisotropic"
  res
}
