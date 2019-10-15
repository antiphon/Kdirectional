#' Kernel Estimate Intensity at Data Points
#'
#' @param x point pattern
#' @param bw bandwidth. Gaussian sd=bw, Epanechnicov domain [-bw, bw].
#' @param kernel Either 'gaussian' or 'epanenchnikov' (partial matching)
#' @param border Border correction to apply.  One of 'none', 'local', 'global', 'toroidal'.
#' @param normalise renormalise so that inverse sum = volume
#' @param loo leave-one-out -estimate? Don't include point i for estimation of int(i)
#'
#' @details
#' For border correction, the bounding box of the coordinates of x will be used, 
#' so at the moment it works only for axis-aligned cuboids.
#'
#' The kernel will be a product of one dimensional kernels,
#' so Epanechnikov in many dimensions will not be truly isotropic but prefers diagonal directions. 
#'
#' @export

intensity_at_points <- function(x, bw, kernel = "gaussian", border = "local", normalise = FALSE, loo = FALSE){
  
  kernel_i <- pmatch(kernel, c("gaussian", "epanechnikov")) - 1
  border_i <- pmatch(border, c("none", "local","global", "toroidal")) - 1
  
  if(! kernel_i %in% 0:1) stop("kernel input failure")
  if(! border_i %in% 0:3) stop("border input failure")
  if(! as.integer(loo) %in% 0:1) stop("loo input failure")
  x <- check_pp(x)
  xy <- x$x
  bbox <- x$bbox
  # go
  out <- int_at_points_c(xy, bbox, bw, kernel = kernel_i, border = border_i, loo = loo)
  #
  if(normalise){
    V <- prod(apply(bbox, 2, diff))
    k <- V/sum(1/out)
    out <- out/k
  }
  #
  # done
  out
}


#' Kernel Estimate Intensity at non-Data locations
#'
#' @param x point pattern
#' @param loc coordinate matrix for non-data locations at which to estimate the intensity.
#' @param bw bandwidth. Gaussian sd=bw, Epanechnicov domain [-bw, bw].
#' @param kernel Either 'gaussian' or 'epanenchnikov'
#' @param border Border correction to apply.  One of 'none', 'local', 'global', 'toroidal'.
#' @param normalise renormalise so that inverse sum = volume
#' @param loo leave-one-out -estimate?
#'
#' @details
#' For border correction, the bounding box of the coordinates of x will be used, 
#' so at the moment it works only for axis-aligned cuboids.
#'
#' The kernel will be a product of one dimensional kernels,
#' so Epanechnikov in many dimensions will not be truly isotropic but prefers diagonal directions. 
#'
#' @export
intensity_somewhere <- function(x, loc, bw, kernel = "gaussian", border = "local"){
  #
  #
  kernel_i <- pmatch(kernel, c("gaussian", "epanechnikov")) - 1
  if(! kernel_i %in% 0:1) stop("kernel input failure")
  border_i <- pmatch(border, c("none", "local", "global", "toroidal")) - 1
  if(! border_i %in% 0:3) stop("border input failure")
  if(missing(loc)) stop("loc should be a coordinate matrix")
  #
  x <- check_pp(x)
  xy <- x$x
  bbox <- x$bbox
  
  if(ncol(loc) != ncol(xy)) stop("dimension mismatch, ncol(loc) != ncol(data)")
  
  # go
  out <- int_at_anywhere_c(xy, loc, bbox, bw, kernel = kernel_i, border = border_i)
  # done
  out
}


#' Find the optimal smoothing for intensity estimation kernel width
#' 
#' @param x point pattern
#' @param bw_vector Vector of bandwidth values to optimize over
#' @param kernel Either 'gaussian' or 'epanenchnikov'
#' @param ... ignored.
#' @details 
#' 
#' Optimize bandwidth $h$ using the Cronie & van Lieshout method [1]. In short, optimise the loss function
#' \deqn{latex}{\sum 1/\hat\lambda_h(x_i) - |W|)^2}
#' with a kernel intensity estimator. See Reference for the paper. 
#' 
#' Note: The optimisation is done without border correction and without leave-one-out, as suggested in [1].
#' 
#' @return list with elements: 
#' 
#' * opt: the optimal bandwidth, given the inputs
#' 
#' * iopt: index of optimal bandwidth in input vector
#' 
#' * loss: the squared loss function (given above)
#' 
#' * statistic: the statistic which root is to be found 
#' 
#' * area: window area |W| used for optimisation.
#' 
#' 
#' 
#' @references 
#' 1. Cronie O, van Lieshout MNM. Bandwidth selection for kernel estimators of the spatial intensity function. 
#' 2016;1-20. Available from: http://arxiv.org/abs/1611.10221
#' 
#' @export
intensity_bandwidth_profile <- function(x, bw_vector, kernel = "gaussian", ...){
  kernel_i <- pmatch(kernel, kerns <- c("gaussian", "epanechnikov")) - 1
  if(! kernel_i %in% 0:1) stop("kernel input failure")
  
  x <- check_pp(x)
  xy <- x$x
  bbox <- x$bbox
  # a reasonable vector?
  if(missing(bw_vector)) {
    bwM <- max(bbox_sideLengths(bbox))
    bw_vector <- exp( seq(log(bwM * 0.01), log(bwM * 2), l = 30) )
  }
  # compute the estimates
  out <- int_bw_c(xy, bw_vector, kernel = kernel_i)
  #
  # The target statistic of which root is the answer
  V <- bbox_volume(x$bbox)
  invs <- colSums(1/out)
  S <- invs - V
  # the loss
  S2 <- S^2
  # the optimal
  iopt <- which.min(S2)
  out <- list(opt = bw_vector[iopt], iopt=iopt, bw = bw_vector, loss=S2, statistic = S, area = V, kernel = kerns[1+kernel_i], dim = ncol(xy))
  class(out) <- c("intensity_bw", is(out))
  out
}

#' Quick plot for intensity bw object
#' 
#' @param x output of \link{intensity_bandwidth_profile}
#' @param ... passed on to plot.default
#' @export
plot.intensity_bw <- function(x, ...) {
  plot(x$bw, x$statistic, xlab="bandwidth", ylab="statistic", ...)
  abline(h=0)
  abline(v=x$opt, lty=2)
}


#' Print for intensity bw object
#' 
#' @param x output of \link{intensity_bandwidth_profile}
#' @param ... ignored
#' @export
print.intensity_bw <- function(x, ...) {
  cat("Kernel smoothing bandwidth estimate:")
  cat("\n* Kernel:", x$kernel)
  cat("\n* Opt. bandwidth:", x$opt)
  cat("\n* Loss at optimum:", x$loss[x$iopt])
  cat("\n* Dimension:", x$dim)
  cat("\n* Window area:", x$area)
  cat("\n")
}
