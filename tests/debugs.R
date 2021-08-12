library(sphere)
library(devtools)
load_all(".")

# # s <- readRDS("../../papers/directed_summaries/code/example_patterns/intensity_1000_small_cut_compressed_0.6.rds")
# # x<-s$x_cox3
# # 
# # fry <- fry_ellipsoids(x, eps=0.1, nangles = 3, nvec = (1:15)*3, verbose=TRUE, r_adjust=1, origin=TRUE)
# 
# #ave <- mean_ellipsoids( el <- fry_ellipsoids(x, border=TRUE, eps=.15, nangles = 2)$el   )
# 
# #ave<-mean_ellipsoids(els<-(fry<-fry_ellipsoids(x, border=TRUE, r_adjust=0.5, eps=0.1, nangles=if(ncol(x$x)==3) 3 else 31, nvec = 1:10))$ellipso )
# 
# #confint(fry$el[[2]])
# 
# 
# load("~/work_off_dropbox/clustered_matern_thomas.rda")
# load("~/work_off_dropbox/clustered_matern_thomas_compressed.rda")
# 
# p <- pp$thomas1
# pc <- ppc$thomas1
# 
# x <- list(x=cbind(p$x,p$y), bbox=with(p$wind, cbind(xrange, yrange)))
# 
# fry <- fry_ellipsoids(x, r_adjust = 0.1, nangles=50, verbose=TRUE)
# 
# print(  sum(is.na(fry$fry))  )


x <- z <- readRDS("debug_profile.rds")
grid <- z$grid
r <- z$r
p <- anisotropy_profile(z, grid, r = r)
