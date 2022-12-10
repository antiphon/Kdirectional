
#' Constructor
#' @export
setClass("fryellipsoids")

#'
#'@export
setClass("confint_fryellipsoids")


############################################################

#' Confidence intervals for Fry ellipsoids
#' 
#' @param x Fitted Fry-ellipsoids object
#' @param fun contrast function, passed on to confint.ellipsoid
#' @param ... Passed to confint.ellipsoid
#' 
#' @details 
#' 
#' See the 'confint.ellipsoid' in 'ellipsoid' for details of ...
#' 
#'
#' @export
confint.fryellipsoids <- function(x, fun, ...) {
  els <- x$ellipsoids
  nvec <- x$n
  # compute the confidence intervals
  ci <- lapply(els, confint, fun=fun, ...)
  out <- list(ci=ci, nvec = nvec, dim = x$dim, own_contrast = !missing(fun))
  class(out) <- "confint_fryellipsoids"
  out
}


#' Summary method for confint_fryellipsoid
#' 
#' @param Fitted Fry-ellipsoids object
#' @param ... ignored
#' 
#' @details 
#' 
#' Extracts the CI-table.
#' 
#' @export
summary.confint_fryellipsoids <- function(x, ...) {
  ci <- x$ci
  r <- nrow(ci[[1]])
  tab <- t(sapply(ci,  function(ta) ta[r,]  ))
  rownames(tab) <- paste0("Contour_", x$nvec)
  tab
}

#' Print method for confint_fryellipsoid
#' 
#' @export
print.confint_fryellipsoids <- function(x, ...) {
  nvec <- x$nvec
  cat("Confidence intervals for", length(nvec), "cumulative count contours estimated from Fry points.\n")
  #
  # get the most interesting one
  tab <- summary(x, ...)
  type <- if(x$own_contrast) "provided contrast statistics"
          else "Equality of semi-axes"
  cat("\n", type, "\n")
  print(tab)
}


#' Plot confidence intervals
#' 
#' 
#' @export
plot.confint_fryellipsoids <- function(x, ylim, ...) {
  tab <- summary(x, ...)
  nvec <- x$nvec
  ci <- tab[,c(4,1,5)]
  if(missing(ylim)){
    ylim <- range(ci)
    ylim <- ylim + c(-1,1)*diff(ylim)*0.05
  }
  #
  plot(NA, xlim=c(0, max(nvec)+1), type="b", 
       xlab="Contour", 
       ylab="", 
       main="", ylim=ylim)
  abline(h=0)
  plot1 <- function(i) { lines(c(nvec[i], nvec[i]), ci[i, c(1,3)], lwd=1) }
  for(i in 1:length(nvec)) plot1(i)
  points(nvec, ci[,2], pch=15)
}



