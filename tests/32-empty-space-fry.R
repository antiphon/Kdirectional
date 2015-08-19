# Testaa tyhjio-fry

library(devtools)
load_all(".")
#load("~/work_off_dropbox/clustered_matern_thomas.rda")
#load("~/work_off_dropbox/clustered_matern_thomas_compressed.rda")

bbox <- with(pp$matern1$window, cbind(xrange, yrange))
bboxc <- with(ppc$matern1$window, cbind(xrange, yrange))

n0 <- 200
nx <- 20
ny <- nx
nxc <- nx*diff(bboxc[,1])/diff(bbox[,1])
nyc <- ny*diff(bboxc[,2])/diff(bbox[,2])
r <- 1
#u <- as.matrix(expand.grid(x=seq(bbox[1,1]+r, bbox[2,1]-r, length=nx),
#                           y=seq(bbox[1,2]+r, bbox[2,2]-r, length=ny)))   
u <- apply(bbox, 2, function(ab) runif(n0, ab[1]+r, ab[2]-r))
#uc <- as.matrix(expand.grid(x=seq(bboxc[1,1]+r, bboxc[2,1]-r, length=nxc),
#                           y=seq(bboxc[1,2]+r, bboxc[2,2]-r, length=nyc)))   #apply(bbox, 2, function
uc <- apply(bboxc, 2, function(ab) runif(n0, ab[1]+r, ab[2]-r))

p <- pp$thomas1
pc <- ppc$thomas1
x <- rbind( cbind(p$x,p$y), u )
xc <- rbind( cbind(pc$x,pc$y), uc )
n <- p$n
nc <- pc$n

N <- n+n0
Nc <- nc + n0
nd <- nrow(u)
ndc <- nrow(uc)

f <- pairwise_vectors(x, from = 1:nd+n, to=1:n, as.xyz = TRUE)
fc <- pairwise_vectors(xc, from = 1:ndc+nc, to=1:nc, as.xyz= TRUE)
#r <- sqrt(sort(rowSums(f^2)))
#rc <- sqrt(sort(rowSums(fc^2)))


par(mfrow=c(3,2))
plot(p, asp=1); points(u, col=2, cex=.4)
plot(pc, asp=1); points(uc, col=2, cex=.4)
cc <- function() points(0,0, pch=4, col=3)
l0 <- c(-1,1) * 8
plot(f, pch=".", asp=1, xlim=l0, ylim=l0);cc()
plot(fc, pch=".", asp=1, xlim=l0, ylim=l0);cc()
l <- c(-1,1)*.25
plot(f, pch=1, cex=.3, asp=1, xlim=l, ylim=l);cc()
plot(fc, pch=1, cex=.3, asp=1, xlim=l, ylim=l);cc()


# not really good