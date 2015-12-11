## ---- include=FALSE------------------------------------------------------
library(spatstat)
library(Kdirectional)
library(ellipsoid)
library(sphere)

## ---- eval=FALSE---------------------------------------------------------
#  library(spatstat)
#  library(Kdirectional)
#  library(ellipsoid)
#  library(sphere)

## ------------------------------------------------------------------------
comp <- 0.7
D <- diag( c(comp, 1/comp))
angle <- pi/8
R <- sphere::rotationMatrix(az=angle)[-3,-3]
C <- R%*%D

## ---- warning=FALSE------------------------------------------------------
set.seed(1)
x0 <- rStrauss(500, gamma = 0.05, R = 0.03, W = affine(square(), solve(C)))

## ------------------------------------------------------------------------
x <- affine(x0, C)

## ---- fig.width=8--------------------------------------------------------
plotp <- function(...) plot(..., cex=.3)
par(mfrow=c(1,2), mar=c(1,1,1,1))
plotp(x0)
plotp(x)

## ---- fig.height=5, fig.width=5------------------------------------------
nna <- nnangle(coords(x))
ang <- nna[,1]-pi
angles.flower(ang, ci = T, k=12)

## ------------------------------------------------------------------------
f <- fry_ellipsoids(x, nvec=1:20, nangles = 30, eps=0.1, r_adjust = .2)

## ---- fig.height=5, fig.width=5------------------------------------------
plot(f, zoom=.3, used_points = FALSE)

## ---- fig.width=5, fig.height=5------------------------------------------
ci <- confint(f, nsim=500)
print(summary(ci)[1:5,], digits=2)
plot(ci, ylim = c(-1,1)*.1)

## ------------------------------------------------------------------------
contours <- f$ellipsoids[-1]
mean_e <- mean_ellipsoids(contours, just_rotation = TRUE)
summary(mean_e)

## ------------------------------------------------------------------------
c(true=angle, est=mean_e$rot_angle + pi)

## ------------------------------------------------------------------------
Rhat <- sphere::rotationMatrix(az=pi)[-3,-3]%*%mean_e$rot 
xr <- affine(x, solve(Rhat))

## ------------------------------------------------------------------------
gamma_vec <- seq(0.5, 1, by=0.05)
grid <- cbind(gamma_vec, 1/gamma_vec)

## ------------------------------------------------------------------------
r <- seq(0, 0.2, length=50)
eps <- pi/10
ani <- anisotropy_profile_fast(xr, grid = grid, r=r, eps=eps, power=2)
ani0 <- anisotropy_profile_fast(x, grid = grid, r=r, eps=eps, power=2)

## ---- fig.width=6, fig.height=6------------------------------------------
plot(ani, type="b", scale=TRUE)
lines(ani0, type="b", col=2)

## ------------------------------------------------------------------------
( s <- summary(ani) )

## ------------------------------------------------------------------------
c(true=comp, est=s$estimate[1])

## ---- fig.width=5, fig.height=3------------------------------------------
Dhat <- s$estimate
xhat <- affine(xr, solve(Dhat))

## ---- fig.width=6, fig.height=5------------------------------------------
par(mfrow=c(2,2), mar=c(1,1,1,1))
plotp(x)
plotp(xr)
plotp(xhat)
plotp(x0)

## ------------------------------------------------------------------------
Chat <- Rhat %*% Dhat
e0 <- ellipsoid_OLS( rellipsoid(1000, c(1,1), R = C) )
e <- ellipsoid_OLS( rellipsoid(1000, c(1,1), R = Chat) )

## ---- fig.height=5, fig.width=5------------------------------------------
plot(e0, add=FALSE)
plot(e, col=3)

