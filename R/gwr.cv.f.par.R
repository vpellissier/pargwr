#' for any given bandwidth, gwr.cv.f.par compute the cv score. Done by running the n GLMs accross multiple cores/nodes
gwr.cv.f.par<-function (bandwidth, y, x, coords, gweight, verbose = TRUE, longlat = FALSE,
                        RMSE = FALSE, weights, show.error.messages = TRUE, ncores)
{
  n <- NROW(x)
  cv <- numeric(n)
  options(show.error.messages = show.error.messages)

  # function used to compute the CV score (leave-one-out approach). For each cell, runs one weigthed GLM per cell, 
  # without including the focal cell, and return observed minus fitted for the focal cell.
  cv.compz<-function(i)
  {
    xx <- x[i, ]
    dxs <- spDistsN1(coords, coords[i, ], longlat = longlat)
    if (!is.finite(dxs[i]))
      dxs[i] <- .Machine$double.xmax/2
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
  
  sfExport(list=c("bandwidth"))
  sfSapply(seq(n), cv.compz)->cv
  
  score <- sum(t(cv) %*% cv) #MSE
  if (RMSE)
    score <- sqrt(score/n) #RMSE
  if (!show.error.messages)
    options(show.error.messages = TRUE)
  if (verbose)
    cat("Bandwidth:", bandwidth, "CV score:", score, "\n")
  score
}