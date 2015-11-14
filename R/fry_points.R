#' Fry points
#' 
#' Pairwise difference vectors of x
#' 
#' @param border Should border effects be cancelled by minus sampling, if TRUE, "double" ignored.
#' @param double Should we return both directions of a pair
#' 
#' @return 
#' 
#' Returns the lengths and direction unit vectors.
#' 
#' @export
fry_points <- function(x, double=FALSE, border=TRUE){
  x <- check_pp(x)
  dim <- ncol(x$x)
  
  if(border){
    n <- nrow(x$x)
    pM <- pairwise_vectors(x$x, asMatrix = TRUE)
    r <- pM$distance
    bdist <- bbox_distance(x$x, x$bbox)
    ok <- do.call(rbind, sapply(1:(n-1), 
                                function(i) {
                                  m<-which(bdist[i] >  r[i,-i]) 
                                  if(length(m))cbind(i,m) else NULL
                                }  
    )  
    )
    # drop diagonal
    ok <- ok[ok[,1]!= ok[,2],]
    # double is implied  
    # convert to xy(z)
    ang <-  if(dim==2) pM$angle[cbind(ok)] else cbind(pM$azi[cbind(ok)], pM$inc[cbind(ok)])
    data_units <- if(dim==2) cbind(cos(ang), sin(ang)) else ai2xyz(ang)
    data_r <- r[cbind(ok)]
    
  }
  else{ # no border correction
    p <- pairwise_vectors(x$x)
    r <- p[,1]
    ang <- p[,-1]
    data_units <- if(dim==2) cbind(cos(ang), sin(ang)) else ai2xyz(ang)
    ##### This is our data for following
    data_r <- r
    data <- data_r * data_units
    if(double){
      data_units <- rbind(data_units, -data_units)
      data_r <- c(data_r, data_r)
      data <- rbind(data, -data)
    }
  }
  list(fry_r=data_r, fry_units=data_units)
}

