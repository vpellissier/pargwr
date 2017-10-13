# Internal function. Computes the local R² for a given point.
# DO NOT RUN OUTSIDE gwr_par()!
local.R2<-function(y, cell, coords, yhat, longlat, adapt, 
	weights, kernel, bandwidth)
{
	i<-cell
    dists.local<-spDistsN1(coords, coords[i,], longlat)
    if(any(!is.finite(dists.local)))
        stop("Infinite distances!")

    if(adapt==TRUE)
        stop("Adaptative bandwidth not implemented yet")

    if(adapt==FALSE){
        if(kernel=="gaussian"){
            weights.i<-weight.gaussian(dists.local, bandwidth)
        }

        if(kernel=="bisquare"){
            weights.i<-weight.bisquare(dists.local, bandwidth)
        }

        if(kernel=="boxcar"){
            weights.i<-weight.boxcar(dists.local, bandwidth)
        }    
    }

    weights.i<-weights.i*weights

    RSS <- sum(weights.i * (y - yhat)^2)
    yss <- sum(weights.i * (y - weighted.mean(y, weights.i))^2)
    return(1 - (RSS/yss))
}