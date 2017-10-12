#' Internal function. Computes the local LM for a given point.
#' 
#' For internal use only
#' 
#' @param x x values
#' @return parameters of a local LM
#' @export
gwr.internal<-function(x, y, cell, coords, bandwidth, weights=NULL,
    kernel, longlat, adapt, se.fit, diagnostic)
{
    i<-cell
    dists.local<-spDistsN1(coords, coords[i,], longlat)
    if(any(!is.finite(dists.local)))
        dists.local[which(!is.finite(dists.local))] <- 0

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

    sum.weights<-sum(weights.i)

    lm.i<-lm.wfit(x, y, w=weights.i)

    coeffs.i<-lm.i$coefficients
    pred.i<-sum(x[i,] * coeffs.i)
    gwr.e.i<-lm.i$residuals[i]

    if(diagnostic==TRUE){
        if(rank(ml.i)!=ncol(x))
            warning("Local LM not full rank")
        invZ<-chol2inv(lm.i$qr$qr[1:rank(lm.i), 1:rank(lm.i)])
        lhat.i<-t(x[i,]) %*% invZ %*% t(x) %*% diag(weights.i)

    }

    if(se.fit==TRUE) 
        coeffs.se.i<-diag(invZ)

    df.i<-c(sum.weights, coeffs.i, coeffs.se.i, pred.i, gwr.e.i)
    names(df.i)<-c("sum.weights", names(coeffs.i), 
                    paste0("SE_", names(coeffs.i)), "yhat", "gwr.error")
    
    return(list(df.i=df.i, lhat.i=lhat.i))
}
