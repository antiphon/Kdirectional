#' Plot a typical conical directed function of range
#' 
#' @param x Output from some function
#' @param r_scale Plot with x-axis r*r_scale
#' @param rmax plot upto this range
#' @param ylim optional range for y-axis
#' @param legpos legend position
#' @param ... passed on to plot
#' 
#' @details Used by: nn_directional_distance_* 
#' 
#' @export

plot.f_cone <- function(x, r_scale=1, rmax, ylim, legpos="topright", lwd = 1, ...) {
  # cut r
  if(!missing(rmax)) x <- x[x$r<rmax,]
  if(missing(ylim)) ylim <- range(x[,-which(names(x)=="r")], na.rm=T)
  if(all(is.na(x$theo))) x$theo <- NULL
  #
  fname <- attr(x, "fname")
  plot(x$r*r_scale, x[,2], col=1, xlab="r", 
       ylab=if(is.null(fname)) "pcf_anin" else fname, type="l", lty=1, ylim=ylim, lwd = lwd, ...)
  n <- ncol(x)
  for(i in 3:n){
    lines(x$r*r_scale, x[,i], col=i-1, lty = i-1, lwd = lwd)
  }
  legend(legpos, names(x)[-1], lty=1:(n-1), col=1:(n-1))
}






