# #' (azi,incl) to 3d coordinates
# #' @export
# ai2xyz <- function(aziinc) {
#   azi <- aziinc[,1]
#   inc <- aziinc[,2]
#   wx<-sin(inc)*cos(azi)
#   wy<-sin(inc)*sin(azi)
#   wz<-cos(inc)
#   cbind(x=wx,y=wy,z=wz)
# }

# In sphere