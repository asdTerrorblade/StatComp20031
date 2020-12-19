#' @title Spearman correlation
#' @description Calculating Spearman correlation
#' @param x the first ordinal data(vector)
#' @param y the second ordinal data(vector)
#' @return value of Spearman correlation
#' @examples
#' \dontrun{
#' x <- runif(10,1,10);y <- runif(10,1,10)
#' x <- ceiling(x); y <- ceiling(y)
#' corS(x, y)
#' }
#' @export
corS <- function(x, y) {
  n <- length(x)
  r <- 1 - 6 * sum((x - y)^2) / (n^3 - n)
  r
}

#' @title Kendall correlation
#' @description Calculating Kendall correlation
#' @param x the first ordinal data(vector)
#' @param y the second ordinal data(vector)
#' @return value of Kendall correlation
#' @examples
#' \dontrun{
#' x <- runif(10,1,10);y <- runif(10,1,10)
#' x <- ceiling(x); y <- ceiling(y)
#' corK(x, y)
#' }
#' @export
corK <- function(x, y)
{
  n <- length(x)
  r <- 0
  for (j in 2:n) {
    for (i in 1:(j-1)) {
      r <- r + sign((x[i] - x[j]) * (y[i] - y[j]))
    }
  }
  2 * r/(n * (n - 1))
}

#' @title Gamma correlation
#' @description calculating Gamma correlation
#' @param d two ordinal random variables distribution
#' @return vealue of Gamma correlation
#' @examples 
#' \dontrun{
#' d <- matrix(ceiling(runif(12,1,12)), 3, 4)
#' corG(d)
#' }
#' @export
corG <- function(d) {
  r <- nrow(d)
  c <- ncol(d)
  x <- y <- 0
  for (i in 1 : (r-1)) {
    for (j in 1 : (c-1)) {
      x <- x + d[i,j] * sum(d[c((i+1) : r),c((j+1) : c)])
    }
  }
  for (i in 1 : (r-1)) {
    for (j in 2:c) {
      y <- y + d[i,j] * sum(d[c((i+1) : r),c(1 : (j-1))])
    }
  }
  (x-y) / (x+y)
}