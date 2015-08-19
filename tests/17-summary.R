#'
library(devtools)
load_all(".")

library(rstrauss)

pp <- rstrauss(1000, bbox=cbind(0:1, 0:1, 0:1))

x <- nnangle(pp$x)
#' drop edge
ok <- bbox_distance(pp$x, pp$bbox) > x[,2]
s <- summary(x)
