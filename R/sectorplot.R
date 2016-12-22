#' Plot a sector estimate
#' 
#' @param x result object from one of the directed summaries
#' @param res smoothness of arc
#' @param col function to derive the colors, see values2colors
#' @param zlim passed to values2colors
#' @param overlapfactor <1 leaves a bit of space between sectors
#' @param ... passed to plot, mostly the "main" parameter
#' 
#' @export

sectorplot <- function(x, res=4, col=rainbow, zlim, overlapfactor=0.95, ...){
  # Ok: 
  ok <- c("pcf_sector", "sector", "Gsector")
  if(!length(intersect(ok, is(x))))stop("Input type not supported by sectorplot.")
  if(x$dim!=2)stop("Only for 2d.")
  # Get the sectors
  ang <- x$theta[[1]]
  eps <- x$epsilon
  sang <- unname( t( rbind(ang, ang) + c(-1,1)*overlapfactor*eps) )
  
  antipode <- "pcf_sector"%in%is(x)
  
  # get the estimates
  val <- x$est
  m <- max(x$r)
  if(missing(zlim)) zlim <- range(x$est, na.rm=T)
  r <- unname( x$r )
  
  plot(NA, xlab="", ylab="", xaxt="n", yaxt="n", asp=1, xlim=c(-m,m), ylim=c(-m,m), ...)
  
  for(i in 1:ncol(val)){  
    co <- values2colors(val[,i], col=col, zlim=zlim)
    draw.arc(0,0, radius=r, angle1 = sang[i, 1], angle2=sang[i, 2], col=co, lwd = 2, n=res)
    if(antipode) draw.arc(0,0, radius=r, angle1 = sang[i, 1]+pi, angle2=sang[i, 2]+pi, col=co, lwd = 2, n=res)
  }
}

