# function used to compute the CV score (leave-one-out approach). For each cell i, runs one weigthed GLM per cell, 
# without including the focal cell, and return observed minus fitted for the focal cell.

cv.compz<-function(i)
{
  xx <- x[i, ]
  dxs <- spDistsN1(coords, coords[i, ], longlat = longlat)
  if (!is.finite(dxs[i]))
    dxs[i] <- .Machine$double.xmax/2
  
  if(kernel=="gaussian")
    w.i<-weight.gaussian(dxs, bandwidth)
  
  if(kernel=="bisquare")
    w.i<-weight.bisquare(dxs, bandwidth)
  
  w.i <- gweight(dxs^2, bandwidth)
  w.i[i] <- 0
  w.i <- w.i * weights
  if (any(w.i < 0 | is.na(w.i)))
    stop(paste("Invalid weights for i:", i))
  lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
  if (!inherits(lm.i, "try-error")) {
    b <- coefficients(lm.i)
    return(weights[i] * y[i] - (t(b) %*% (weights[i] *
                                            xx)))
  }
}
