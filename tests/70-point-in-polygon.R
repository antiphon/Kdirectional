# Point inside polygon, remove dependency on sp?

# NOTE: Not removed, code in sp is quite non-trivial

library(devtools)
load_all()

library(rgl)
# The code
triangle <- rbind( c(0, 0, 0), c(1, 1, 0), c(0, 1, 1) )
pl <- three_points_to_plane(triangle) |> cbind()
d <- 3

n <- ncol(pl)
bi <- 1:d
ai <- bi+d
linesf <- list(lines, if(exists("lines3d"))lines3d else NULL) [[d-1]]
for(i in 1:n) linesf((rbind(pl[ai,i], pl[bi,i]+pl[ai,i])), col = 2)






# Random points in sphere
u     <- matrix(runif(30, -1, 1), ncol = 3)
u     <- u/apply(u, 1, \(u) sqrt(sum(u^2)))

proj  <- project_to_plane(u, pl[1:3])
tproj <- project_to_plane(triangle, pl[1:3])
# check
cbind(sp  = ii <- sp::point.in.polygon(proj[,1], proj[,2], tproj[,1], tproj[,2]),
   my = point.in.polygon(proj[,1], proj[,2], tproj[,1], tproj[,2]) )



#####
plot3d(triangle, aspect = FALSE, type = "s", size = 4, col = 1:3)
triangles3d(triangle, alpha = .2)
spheres3d(u, col = c("yellow", "blue")[1+ii], radius = 0.05)
