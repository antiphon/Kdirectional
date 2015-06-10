# #' values to colors
# #' 
# #' Turn a vector of values to colors.
# #' 
# #' 
# #' 
# #' @export
# values2colors <- function(v, n=100, zlim, col=heat.colors, ...){
#   if(missing(zlim)) zlim <- range(v)
#   pr <- v
#   pp <- (pr-zlim[1])/diff(zlim)
#   e <- floor(pp*(n-1)) + 1
#   hc <-col(n)
#   hc[e]
# }

# Sphere