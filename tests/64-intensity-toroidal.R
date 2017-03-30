# Int with toroidal distances

library(devtools)
load_all(".")


set.seed(1)
# set dimension:
dim <- 2

# bbox
lx <- c(0,3.1)
ly <- c(-.2, 2)
bbox <- cbind(lx, ly, lx)[,1:dim]
V <- prod(apply(bbox, 2, diff))

# number of points
N <- 200
int <- N/V

# sample unif pattern
x <- apply(bbox, 2, function(ab) runif(N, ab[1], ab[2]))
p <- list(x=x, bbox = bbox)


# bandwidth
bw <- .15

# optimize bw?
if(0){
  area <- function(x) V
  coords <- function(x)x
  opt <- intensity_bandwidth_profile(x, seq(0.05, 5, l = 30))
  bw <- opt$opt
  par(mfrow=c(1,2))
  plot(opt$loss, log="y")
}

# Estimate lambda
# at grid
ng <- if(dim==2) 100 else 50
sx <- seq(bbox[1,1], bbox[2,1], l = ng)
sy <- seq(bbox[1,2], bbox[2,2], l = ng)
if(dim==3) sz <- seq(bbox[1,3], bbox[2,3], l = ng)
grid <- if(dim==2) as.matrix( expand.grid(sx, sy) ) else as.matrix( expand.grid(sx, sy, sz) )

g <- intensity_somewhere(p, loc = grid, bw = bw, border = bo <- "t")

# at points
v <- intensity_at_points(p, bw = bw, border = bo)

# nearest pixel
ni <- apply(x, 1, function(v)  which.min(colSums((t(grid)-v)^2)  ))

# Plot
G <- if(dim==3) array(g, dim=c(ng,ng,ng))[,,ng/2] else matrix(g, nrow=length(sx))
xinc <- if(dim == 3)  which(abs(x[,3]-sz[ng/2]) < 0.1) else 1:N


co <- heat.colors(1200)
zl <- range(G, v[xinc])
vz <- (v-zl[1])/diff(zl)
cox <- co[ 1 + round(1199 * vz) ]


par(mfrow=c(2, dim))
image(sx, sy, G, asp=1, col = co, zlim = zl)
points(x[xinc,1:2], col= cox[xinc], pch=19)
points(x[xinc,1:2], pch=".")

plot(v-g[ni], col=4, pch=19, cex=.5, main="difference lambda(x) - lambda(grid(x))")
abline(h=0)
for(i in 1:dim)
  plot(x[,i], v, main=paste0("along dim ", i), ylab="intensity(x)")

