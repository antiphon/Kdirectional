#' plotting of 3d ellipses
library(devtools)
load_all(".")

x <- rellipsoid(1000)
e <- ellipsoid_OLS(x)
plot(e)
