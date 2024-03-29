# Dev and Test the rosenberg wavelet transform
library(devtools)
load_all(".")
Rot <- Kdirectional::rotationMatrix3(az=pi/2)[-3,-3]
#source("tests/test-data.R")

#x <- pp

nn <- 100
x <-  cbind(runif(nn), runif(nn))  * 60
#x[,2] <- round(x[,2])



#source("~/Sync/Aka/julkaisut/anisotropy-review-overleaf/examples/synthetic-examples/generate-stripe.R")
#x <- stripes


x <- check_pp(x)
bb <- x$bbox

res <- 180

# End results
if(1){

  b <- rosenberg(x, steps = res, include = 0.8, theta = seq(0,pi, l = 180))
  b0 <- rosenberg(x, steps = res, include = 1, theta = b[,1])
  ##
  par(mfrow = c(2,2))
  plot(x$x, asp = 1)
  plot(b, type="l", xaxt="n")
  lines(b0, col = 2)
  axis(1, at = seq(0, pi, l = 5), labels = expression(0, pi/4, pi/2, 3*pi/4, pi))
  plot(pcf_anisotropic(x))
}

## Check the counts per direction, something wrong with along-the-x-axis
if(0){
  co <- c_rosenberg_intensities(x$x, x$bbox, res)
    
}


# envelope?

if(0){
  library(spatstat)
  asfv <- function(x, ...) as.fv(rosenberg(x, steps = res))
  z <- envelope( ppp(x$x[,1],x$x[,2], window = as.owin(c(x$bbox))) , asfv)
  plot(z)  
}

## check the normals
## ok! obsolete
if(0){
  print(b)
  print( bbquad_planes(bbox2bbquad(bb)  ))
  #print( unclass(bbox2bbquad(bb)))
}


## line hit planes 
if(0){
  bq <- bbox2bbquad(bb)  
  planes <- bbquad_planes(bq)
  a <- 2.199115  #a <- pi/2#runif(1, 0, pi)
  v <- 14
  line   <- c(x1 = x$x[v,1], x2 = x$x[v,2], u1 = cos(a), u2 = sin(a))
  hit <- line_hit_planes(line, planes)
  plot(bq$vertices, pch=NA, asp = 1 , xlim = c(-1,2), ylim = c(-1,2))
  plot.bbquad(bq, add = TRUE)
  points(rbind(line[1:2]))
  lines(rbind(line[1:2] - 3*line[3:4], line[1:2] + 3*line[3:4]))
  points(t(hit), pch= 4)
}



# check counts and "sector area inside box" approximation
#OK  OBSOLETE
if(0){
  i <- 10
  z<-x$x[i,,drop=FALSE]
  b <- c_rosenberg_intensities(x$x, x$bbox, res)
  plot.bbquad(bbox2bbquad(bb), add = FALSE, asp = 1)
  for(ai in 0:(res-1) ){
    d <- 1/res * pi 
    ang <- ai/res * pi #- d/2
    # sectors of size 1deg
    p0 <- cbind(cos(ang), sin(ang))*100
    p1 <- cbind(cos(ang+d), sin(ang+d)) * 100
    pol <- rbind(z, p0, p1, z, -p0, -p1)
    #col <- b[i,ai+1]+1  
    col <- rgb(0,0,b[i,ai+1] / max(b[i,]))
    polygon(pol , col = col, 
            border=NA )
  }
  plot.bbquad(bbox2bbquad(bb), add = TRUE, asp = 1, ecol =  "white")
  points(z, pch=19, col = "yellow")
  points(x$x, pch=3, cex = .8, col ="white")
  
}


