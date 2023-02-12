#' Fit second order intensity ellipses or ellipsoids
#' 
#' @param x point pattern
#' @param nvec which jumps to consider
#' @param r_adjust consider only fry points with length r_adjust*max(window side) 
#' @param eps add to the sector angle. Results in smoothing if >0 (default=0).
#' @param nangles If 3d, number of refinements (def=2), in 2d def=20 sectors.
#' @param cylindric Use a cylinder instead of a sector
#' @param border Should border effects be cancelled by minus sampling, if TRUE, "double" ignored.
#' @param double Should we use both directions of a pair
#' @param origin Constrain the ellipses' centers to origo.
#' @param keep_data Keep the data for each ellipsoid, inside each ellipsoids? Helps with averaging.
#' @param data The vectors to which fit ellipses. Default is missing, and we compute Fry-points.
#' @param ALS (false) Use ALS to refine OLS fits.
#' 
#' @details 
#' 
#' First we compute the Fry-points, then we look at the Fry-points in many directions from 
#' origin and find the kth point in each, then we fit ellipsoids to the kth-points.
#' 
#' @export

fry_ellipsoids <- function(x, nvec=1:5, r_adjust=1, nangles, eps=0, 
                            cylindric=FALSE, double=FALSE, 
                            border=TRUE, origin=TRUE, verbose=FALSE,
                            keep_data=TRUE, double2=FALSE,
                            data, ALS=FALSE, ...) {
  x <- check_pp(x)
  dim <- ncol(x$x)
  #
  cat2 <- if(verbose) cat else function(...)NULL
  #
  ## Border correction
  if(missing(data)){
    fry <- fry_points(x, double=double, border=border)
    data_r <- fry$fry_r
    data_units <- fry$fry_units
  }
  else{ # data given
    data_r <- sqrt(rowSums(data^2))
    data_units <- data/data_r
  }
  # drop long range for easier computations
  r_ok <- data_r < r_adjust * max(bbox_sideLengths(x$bbox))
  data_r <- data_r[r_ok]
  data_units <- data_units[r_ok,]
  data <- data_units * data_r
  
  cat2("Fry points gathered, n=", nrow(data),"\n")
  #
  ##### The directions -grid
  if(missing(nangles)) {
    n <- nrow(x$x)
    # use poisson approximation and requiremnt P(fry's in sector) ~= n(x)/const
    nangles <- if(dim==2) round( n/6 )
               else round( n * pi / 20 )
    if(dim == 3) nangles <- (1:5)[ which.min(abs(nangles- 20*4^(1:5)/2+2) )]
  }
  if(dim==2){
    angles <- seq(0, 2*pi, length=nangles+1)[-1]
    delta <- diff(angles[1:2])/2
    angles <- angles -  delta
    grid_unit <- cbind(cos(angles), sin(angles))
  }
  else{
    grid_unit <- directions <- ll2xyz( sphere.grid(N=nangles, ico=TRUE) )
    nng <-  apply(as.matrix(dist(directions)), 1, function(v) sort(v)[2])
    mdi <- mean(nng) # mean distance between directions in the grid
    # angle between two unit directions with distance mdi
    delta <- 2*asin(mdi/(sqrt(3)) )
  }
  if(cylindric){
    # go for the width of the cylinder grown up to nearest neighbour of origin
    # in any direction
    m <- sort(data_r)[1]
    delta <- sin(delta) * m 
  }
  cat2("Direction grid generated, n=", nrow(grid_unit),"\n")
  #
  #### sector width. eps >0 results in overlaps
  delta <- eps + delta
  #
  #### Compute the inclusions per fry point in each cone:
  # conical
  if(!cylindric){
    inside <- apply(grid_unit, 1, 
                    function(u) {
                      s <- data_units%*%u
                      a <- acos(s)
                      a <- pmin(pi-a, a)
                      if(double & !double2) a < delta else a< delta & s>0
                    })
  }
  else { # cylindrical
    inside <- apply(grid_unit, 1, 
                    function(u) {
                      s <- c(data%*%u)
                      p <- data - t(sapply(s, "*", u))
                      d <- sqrt(rowSums(p^2))
                      if(double) d < delta & s >= 0
                      else d < delta
                    })
  }
    ## get the ranges of each fry point, sorted
  Kd     <- apply(inside, 2, function(i) cbind(r=sort(data_r[i])) )
  #
  ns <- sapply(Kd, length) # counts per sector
  #
  cat2("Sector inclusions resolved, range of point count varies ", 
       min(ns), "-", max(ns),"\n")
  # warnings in case data too sparse
  if(max(nvec) > min(ns)) 
    cat2("Requested neighbour count not available for all directions (min is", 
                 min(ns),")\n")
  #
  #
  ## Compute directional isocurves, the ranges of jumps in each direction
  iso_ranges <- sapply(Kd, function(k) {
      ok <- nvec[nvec <= length(k)]
      e <- k[ok, 1]  # need nvec out of this
      n <- length(e)
      if(n < length(nvec)) e <- c(e, rep(NA, length(nvec)-n))
      e
    } )
  #
  #
  # Fit the ellipses:
  iso_ranges <- rbind(iso_ranges)
  # check counts
  counts <- apply(iso_ranges, 2, function(v) sum(!is.na(v)) )
  cat2("Range vectors per direction jumps computed, count varies ", 
       min(counts), "-",max(counts),"\n")
  perc <- apply(iso_ranges, 1, function(v)sum(!is.na(v)))
  cat2("Data points per contour varies ", min(perc), "-",max(perc),"\n")
  #
  good_k <- perc > (3+dim)
  if(any(perc< (3+dim))) 
    cat2(sum(!good_k),"contours not computed, not enough data.\n")
  #
  # fitting of ellipses:
  els <- apply(iso_ranges[good_k, ], 1, function(r){
    p <- r * grid_unit
    p <- p[!is.na(r),]
    if(double2) p <- rbind(p, -p) # in case add antipodes
    el <- ellipsoid_OLS(p, origin=origin) # in 'ellipsoid' package 
    # refine with ALS?
    if(ALS) el <- ellipsoid_ALS(p, s2 = seq(el$ols_fit$s2*.1, el$ols_fit$s2, length=10))
    if(keep_data) el$data <- p
    el
  } )
  #
  nvec <- nvec[good_k]
  # done
  res <- list(ellipsoids=els, dim=dim, param=list(cylindric=cylindric, delta=delta, eps=eps),
              jumps=iso_ranges, 
              fry=data,
              n=nvec, grid_unit=grid_unit, double2= double2, double=double)
  #
  class(res) <- "fryellipsoids"
  res
}

#############################################################

#' plot method for fryellipsoids
#' 
#' @export
plot.fryellipsoids <- function(x, ellipsoids=TRUE, used_points=TRUE, sectors=NULL, xlim=NULL, ylim=NULL, zoom = NULL, pch = 1, cex=0.8, ...) {
  
  
  # zoom
  if(is.null(xlim)) xlim <- range(x$fry[,1])
  if(is.null(ylim)) ylim <- range(x$fry[,2])
  if(!is.null(zoom)) {
    xlim <- zoom * xlim
    ylim <- zoom * ylim
  }
  

  if(x$dim==2){
    plot(x$fry, asp=1, xlab="x", ylab="y", pch=pch, 
         main=paste0("Fry points", ifelse(used_points, ", estim. points(x) "," "), "and fitted ellipses"), 
         cex=cex, xlim=xlim, ylim=ylim)
    for(i in 1:length(x$ellipsoids)){
      if(used_points){
        if(is.null(x$v2)){
          points(x$grid_unit * x$jumps[i,], pch=20, col=i)
        }
        else if(!is.null(x$ellipsoids[[i]]$data)) {
          points(x$ellipsoids[[i]]$data, pch=20, col=i)
        }
      }
      if(ellipsoids) plot(x$ellipsoids[[i]], col=i, ...)
    }
    abline(h=0, col="gray50")
    abline(v=0, col="gray50")
    
    # 
    if(any(sectors)){
      delta <- x$param$delta
      r <- max(x$fry)
      for(i in sectors){
      u <- x$grid_unit[i,]
      if(x$param$cylindric){
        ut <- u[2:1]*c(-1,1)
        uv <- sapply(seq(0, r, length=2), "*", u )
        lines(rbind(delta*ut, -delta*ut))
        lines(t(delta*ut+uv))
        lines(t(-delta*ut+uv))  
      }
      else{
        m <- cbind(c(cos(delta), sin(delta)),c(-sin(delta), cos(delta)))
        u1 <- c(m%*%u)
        u2 <- c(solve(m)%*%u)
        lines(rbind(c(0,0), u1))
        lines(rbind(c(0,0), u2))
      }
      }
    }
    
  }
  else warning("Plot not implemented in 3d")
}

#############################################################
#' Print method for fryellipsoids
#' 
#' @export
print.fryellipsoids <- function(x, ...) {
  cat("Ellipsoid fit to fry points.\n")
  cat("Number of contours:", length(x$n), "\n")
  cat("Directions:", nrow(x$grid_unit), "\n")
  if(x$param$eps) cat("Direction smoothing:", x$param$eps, "\n")
}







