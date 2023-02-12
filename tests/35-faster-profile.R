## Faster profiling: Effect of RC on translation weights.

library(devtools)
load_all(".")

load("test_strauss_3d.rda")
l <- c(-1,1)*.5

# start
pp <- list(x=x, bbox=cbind(l,l,l))
# compress
co <- 0.8; s <- 1/sqrt(co)
M <- diag(c(s,s,co))

y <- x%*%t(M)

bby <- pp$bbox%*%t(M)
ppy <- list(x=y, bbox=bby)
# translation weights

ww <- translation_weights(pp)
wwy <- translation_weights(ppy)

# any change
all.equal(ww,wwy)

# ok, it is the volume presering compression.

# how about border distance
bd <- bbox_distance(pp$x, pp$bbo)
bdy <- bbox_distance(ppy$x, ppy$bbo)
#
# there might not be any faster way.