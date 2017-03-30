# # Anisotropy profile assuming compression in one dimension
# # 
# # Given a direction i, compute the anisotropy profile for set a compressions.
# # 
# # @param x point pattern, list(x=coords, bbox=bounding-box)
# # @param i direction of compression.
# # @param cvec vector of compressions to go through
# # @param ... passed on to anistropy.abs
# # @export
# # 
# 
# anisotropy_profile <- function(x, i=1, cvec, verb=TRUE, ...) {
#   stop("obsolete function, use anisotropy_profile_fast")
#   profile <- NULL
#   cat2 <- if(verb) cat else function(...)NULL
#   if(is.null(x$bbox)) stop("x should be list(x=coords, bbox=bounding-box)")
#   dim <- ncol(x$x)
#   dxyz <- rep(1, dim)
#   # inverse transform the pattern
#   f <- function(ce) {
#     dxyz[i] <- 1/ce
#     dxyz[-i] <- ce^(1/(dim-1))
#     A <- diag(dxyz)[1:dim, 1:dim]
#     bnod <- as.matrix(expand.grid(as.data.frame(x$bbox)))%*%A
#     list(x=t(t(as.matrix(x$x) %*% A)), bbox=apply(bnod,2,range))
#   }
#   k<-0
#   vals <- list()
#   for(ce in cvec) {
#     # deform
#     xd <- f(ce)
#     vals[[paste0("c",ce)]] <- anisotropy.abs(xd, ...)
#     profile <- c(profile, vals[[paste0("c",ce)]]$statistic)
#     cat2(k<-k+1, "/", length(cvec), "   \r")
#   }
#   cat2("\nDone.")
#   profile <- data.frame(compression=cvec, anisotropy=profile)
#   
#   list(profile=profile, raw=vals)
# }
