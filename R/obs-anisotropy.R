# #' Anisotropy metric, absolute value
# #' 
# #' From Redenbach et al 2009
# #' Idea: the data is squeezed in one dimension (expanded in others).
# #' Find this out by discrepancy w.r.t other direction(s).
# #' 
# #' @param x point pattern
# #' @param target If given, compare the average of other directions  to this. i=1,2,3 as x,y,z
# #' @param ... passed on to Kest_directional, will use parameter u.
# #' 
# #' @export
# 
# anisotropy.abs <- function(x, estimates, fun="K", target=NULL, ...){
#   funs <- list(K=Kest_along_axis, F=Fest_along_axis)
#   if(missing(estimates)){ # compute the axis-directed summary functions
#     estimates <- funs[[fun]](x, ...)
#   }
#   vals <- estimates$vals
#   r <- estimates$r
#   dim <- estimates$dim
#   #'
#   #' The deviate:
#   dx <- diff(r[1:2])
#   S <- 0
#   if(is.null(target)){
#     for(i in 1:(dim-1))
#       for(j in (i+1):dim) S <- S + sum(abs(vals[i,]-vals[j,]))
#   }
#   else{
#     S <- sum(abs(vals[target,] - colMeans(rbind(vals[-target,]))))
#   }
#   S <- S*dx
#   #' and we are done.
#   list(statistic=S, estimates=vals, r=r, type="integrated absolute differences", func=fun)
# }
# 
# 
# 
# 

# Obsoloete slow 
