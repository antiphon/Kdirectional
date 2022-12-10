## ---- include=FALSE-----------------------------------------------------------
library(spatstat)
library(Kdirectional)

## ---- eval=FALSE--------------------------------------------------------------
#  library(spatstat)
#  library(Kdirectional)

## -----------------------------------------------------------------------------
comp <- 0.7
D <- diag( c(comp, 1/comp))
angle <- pi/8
R <- rotationMatrix3(az=angle)[-3,-3] # internal, not rgl one
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

## ---- fig.height=5, fig.width=5-----------------------------------------------
nna <- nnangle(coords(x))
ang <- nna[,1]-pi
angles.flower(ang, ci = T, k=12)

## -----------------------------------------------------------------------------
f <- fry_ellipsoids(x, nvec=1:20, nangles = 30, eps=0.1, r_adjust = .2)

