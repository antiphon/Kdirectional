# *-anin with new itnensity estimator
library(devtools)
load_all(".")

if(!exists("z"))load("test_strauss.rda")

x <- z

# int:
bw <- 0.2
l <- intensity_at_points(x, bw = bw)

# K:
a <- Kest_anin(x, lambda = l)
b <- Kest_anin(x, lambda_h = bw, kernel = "g")
all.equal(a,b)

#
a <- Kest_anin_cylinder(x, lambda = l)
b <- Kest_anin_cylinder(x, lambda_h = bw, kernel = "g")

all.equal(a,b)

# g
a <- pcf_anin(x, lambda = l)
b <- pcf_anin(x, lambda_h = bw, kernel = "g")
all.equal(a,b)

a <- pcf_anin(x, lambda = l)
b <- pcf_anin(x, lambda_h = bw, kernel = "g")
all.equal(a,b)

