#' Mugglestone and Renshaw's polar spectrums
#'
#' @param x point pattern, list of $x and $bbox
#' @param r range frequencies
#' @param std Standardize coordinates?
#' @param ... ignored
#' 
#' (r,angle)=0 value is set to NA, and returned as $zerov in the return object.
#' 
#' @export

spectrum <- function(x, std=TRUE, scale=TRUE, ...) {
  if(is.null(x$bbox)) stop("x should be list(x=coordinates-matrix, bbox=bounding-box)")
  bbox <- x$bbox
  dim <- ncol(bbox)
  if(dim!=2) stop("Only for 2D now.")
  sidelengths <- apply(bbox, 2, diff)
  n <- nrow(x$x)
  # Create the grid
  r <- 1:30
  theta <- seq(0, 170, by=10)
  polarv <- as.matrix(  expand.grid(r, theta)  )
  pgrid <- t( apply(polarv, 1, function(x) x[1]*c(cos(x[2]*pi/180), sin(x[2]*pi/180))  ) )
  #
  P <- periodogram(x, std=std, pgrid=pgrid, ...)
  # 
  v <- P$v
  Asp <- Rsp <- NULL
  # R-spectrum:
  for(ri in r){
    ok <- ri-1 < polarv[,1]  & polarv[,1] <= ri
    s <- mean(v[ok], na.rm=TRUE)
    Rsp <- c(Rsp, s)
  }
  # angle spectrum:
  for(a in theta){
    ok <- a-5 < polarv[,2]  & polarv[,2] <= a +5
    s <- mean(v[ok], na.rm=TRUE)
    Asp <- c(Asp, s)
  } 
  # Scale to CSR
  if(scale){
    Rsp <- Rsp/2
    Asp <- Asp/2
  }
  #
  list(Rspectrum=data.frame(r=r, v=Rsp),
       Aspectrum=data.frame(theta=theta, v=Asp),
       P=P, std=std)
}