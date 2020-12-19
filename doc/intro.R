## ----eval=FALSE---------------------------------------------------------------
#  corS <- function(x, y)
#  {
#    n <- length(x)
#    r <- 1 - 6 * sum((x - y)^2)/(n^3 - n)
#    return r
#  }

## ----eval=FALSE---------------------------------------------------------------
#  corK <- function(x, y)
#  {
#    n <- length(x)
#    r <- 0
#    for (j in 2:n) {
#      for (i in 1:(j-1)) {
#        r <- r + sgn((x[i] - x[j]) * (y[i] - y[j]))
#      }
#    }
#    return 2 * r/(n * (n - 1))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  corG <- function(d) {
#    r <- nrow(d)
#    c <- ncol(d)
#    x <- y <- 0
#    for (i in 1 : (r-1)) {
#      for (j in 1 : (c-1)) {
#        x <- x + d[i,j] * sum(d[c((i+1) : r),c((j+1) : c)])
#      }
#    }
#    for (i in 1 : (r-1)) {
#      for (j in 2:c) {
#        y <- y + d[i,j] * sum(d[c((i+1) : r),c(1 : (j-1))])
#      }
#    }
#    (x-y) / (x+y)
#  }

