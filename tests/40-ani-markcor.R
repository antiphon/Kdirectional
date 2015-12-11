#' ani mark correlation
#' 

library(devtools)
load_all(".")

load("test_strauss.rda")

x <- z
marks <- exp(rnorm(nrow(z$x), x$x, 0.1))

# check by pcf
u <- rbind(c(0,1))

res <- markcorr_anisotropic(x, marks, directions = u, bw=c(0.01, .5))

#gr <- pcf_anisotropic(x, res$r, theta=list(0), h=res$bw)
#plot(gr$r, gr$est, type="l", col=2)
#with(res, lines(r, pcf, "l"))

with(res, plot(r, mcor, "l", ylim=c(0,2)))

abline(h=1)


f <- fry_points(x)
#################################################################################
# the kernel, say (r=0.1, (0,1)):
# epa <- function(d, bw) {k<-0*d;k[d<bw] <- 0.75 * (1-(d[d<bw]/bw)^2) /bw;k}
# xs <- seq(-.05,.2, length=200)
# grid <- as.matrix( expand.grid(xs,xs) )
# lgrid <- apply(grid, 1, function(v) sqrt(sum(v^2)))
# u1 <- c(1,1)/sqrt(2)
# abw <- res$bw[2]
# ua <- c(cos(pi/4+abw), sin(pi/4+abw))
# ub <- c(cos(pi/4-abw), sin(pi/4-abw))
# ang <- acos((grid%*%u1)/lgrid)
# ang <- abs(ang) #pmin(ang, pi-ang)
# 
# M <- 0
# for(r1 in seq(0.0, 0.15, by=0.025)){
#   dr <- abs(lgrid-r1)#apply(grid, 1, function(v) sqrt(sum((r1*u1-v)^2) ))
#   kr <- epa(dr, res$bw[1])
#   ka <- epa(ang, abw)
#   M <- M + matrix(kr*ka, ncol=length(xs))
# }
# image(M, asp=1, x=xs, y=xs, useRaster=TRUE)
# lines(rbind(c(0,0),ua,u1,ub, c(0,0)))
# points(f$fry_r * f$fry_units, pch=19, cex=0.1)


