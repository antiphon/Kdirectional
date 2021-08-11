#' Plot for pcf_directions
#' 
#' @export
plot.pcf_directions <- function(x, ymax=2,  ...) {
  # 2D
  if(x$dim==2){
    theta <- round(x$r_phi[,2], 6)
    thetas <- unique(round(sort(theta), 6))
    pcf <- split(x$est, theta)
    nl <- length(pcf)
    r <- split(x$r_phi[,1], theta)
    rangs <- sapply(r, range)
    
    plot(NA, xlim=range(rangs), ylim=c(0, ymax), main="pcf in different directions", ylab="")
    for(i in 1:nl){
      lines(r[[i]], pcf[[i]], col=i, lty=ifelse(nl<10, i, 1))
    }
    if(nl<10) legend("topright", paste("angle", thetas), col=1:length(pcf), lty=1:nl)
  }
  else{#' 3D
    
  }
  warning("obsolete function.")
}

#' Plot a angluar plot of the pcf_directions, 2D
#' 
#' @export
flower.pcf_directions <- function(x, scale_points=TRUE, col=rainbow, zlim, ...) {
  co <- values2colors(x$est, col=col, zlim=zlim, ...)
  ce <- if(scale_points) sqrt(x$r_phi[,1]/max(x$r_phi[,1])+.1) else 1
  L <- range(x$directions[,1])
  plot(x$directions, col=co, pch=19, ylim=L, xlim=L, asp=1, cex=ce, ...)
  points(-x$directions, col=co, pch=19, cex=ce)
}

