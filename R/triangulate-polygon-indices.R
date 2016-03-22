#' Triangulate a polygon given by indices
#' 
#' Chooses the first vertex of the first edge as apex
#' 
#' @export
triangulate_indices_cohen_hickey <- function(edges){
  nt <- length(unique( c(edges) ))-2
  if(nt<1) return(NULL)
  if(nt==1) list(edges)
  i <- edges[1,1]
  done <- 1
  j <- edges[1,2]
  idx <- NULL
  for(l in 1:nt){
    row <- setdiff(which( apply(apply(edges, 2, "==", j), 1, any)), done)
    k <- setdiff(edges[row,], j)
    idx <- rbind(idx, c(i, j, k))
    done <- c(done, row)
    j <- k
  }
  idx
}