# Internal function. Computes the predicted values.
# For internal use only
gwr.pred.internal<-function(x, y, coords, cell, newcoords, 
    newdata, bandwidth, weights, yhat,
    kernel, longlat, adapt, se.fit)
{
    x.vars<-colnames(newdata)
    i<-cell
    dists.local<-spDistsN1(coords, newcoords[i,], longlat)
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
    pred.i<-sum(newdata[i,] * coeffs.i)
    #gwr.e.i<-lm.i$residuals[i]

    RSS<-sum(weights.i * (y - yhat)^2)
    yss<-sum(weights.i * (y - weighted.mean(y, weights.i))^2)
    localR2.i<-1 - (RSS/yss)

    if(se.fit==TRUE){
        stop("SE computation for the prediction are not implemented yet!")
    #    coeffs.se.i<-diag(invZ)
    #    df.i<-c(sum.weights, coeffs.i, coeffs.se.i, pred.i, gwr.e.i)
    #    names(df.i)<-c("sum.weights", x.vars, 
    #                paste0("SE_", x.vars), "prediction", "local.R2")
    }
    else{
        df.i<-c(sum.weights, coeffs.i, pred.i, localR2.i)
        names(df.i)<-c("sum.weights", x.vars, "prediction", "local.R2")
        lhat.i<-NA
    }
    
    return(list(df.i=df.i, lhat.i=lhat.i))
}
