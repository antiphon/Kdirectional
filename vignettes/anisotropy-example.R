## ---- include=FALSE-----------------------------------------------------------
library(spatstat)
library(Kdirectional)
library(ellipsoid)
library(sphere)

## ---- eval=FALSE--------------------------------------------------------------
#  library(spatstat)
#  library(Kdirectional)
#  library(ellipsoid)
#  library(sphere)

## -----------------------------------------------------------------------------
comp <- 0.7
D <- diag( c(comp, 1/comp))
angle <- pi/8
R <- sphere::rotationMatrix(az=angle)[-3,-3]
C <- R%*%D

## ---- warning=FALSE-----------------------------------------------------------
set.seed(1)
x0 <- rStrauss(500, gamma = 0.05, R = 0.03, W = affine(square(), solve(C)))

## -----------------------------------------------------------------------------
x <- affine(x0, C)

## ---- fig.width=8-------------------------------------------------------------
plotp <- function(...) plot(..., cex=.3)
par(mfrow=c(1,2), mar=c(1,1,1,1))
plotp(x0)
plotp(x)

