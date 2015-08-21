#' fry ell confindence intervals
#' 
#' #
library(devtools)
library(sphere)
load_all(".")

if(!exists("pp"))source("tests/test-data.R")
ff <- fry_ellipsoids(pp, nvec=c(1,5,10,25))

x <- confint(ff)

summary(x)

print(x)

plot(x)