#' Parse input for marks
#'
#' @param x Input object
#' @param marks potential marks in case x has none.
#' @export

parse_marks <- function(x, marks) {
  
  dunno <- paste("Unable to parse marks from input pattern of type:", is(x)[1])
  m <- NULL
  if(is(x, "ppp")){ # if we have spatstat object
    m <- x$marks
  }
  if(is.null(m) & !is.null(marks)) {
    m <- marks
  }
  if(is.null(m)) stop(dunno)
  m
}
