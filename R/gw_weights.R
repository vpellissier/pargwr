#' Bisquare weights function
#' 
#' Computes a vector of weights using a bisquare kernel
#'  
#' @param dist distance vector between observations and fit point
#' @param bandwidth cutoff distance (i.e. distance after which weight=0)
#' @export
weight.bisquare<-function (dist, bandwidth)
{
  max.d2 <- bandwidth^2
  dist2<-dist^2
  weights <- ifelse(dist2 > max.d2, 0, (1 - (dist2/max.d2))^2)
  return(weights)
}

#' Gaussian weights function
#' 
#' Computes a vector of weights using a bisquare kernel
#'  
#' @param dist distance vector between observations and fit point
#' @param bandwidth bandwidth
#' @export
weight.gaussian<-function(dist, bandwidth)
{
  weights<-exp(-0.5*((dist/h)^2))
}