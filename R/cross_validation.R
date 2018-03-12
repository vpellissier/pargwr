# function used to compute the CV score (leave-one-out approach). For each cell, runs one weigthed GLM per cell, 
# without including the focal cell, and return observed minus fitted for the focal cell.

cv.compz<-function(i, descr, longlat, kernel, bandwidth)
{
    if(is(descr, "big.matrix.descriptor")){
        ndf<-attr(descr, "description")$ncol
        par.df<-bigmemory::attach.big.matrix(descr)
    }
    else{
        ndf<-ncol(descr)
        par.df<-descr
    }
    
    posy<-c(1)
    posx<-c(2:(ndf-3))
    posweights<-c(ndf-c(2))
    poscoords<-c(ndf-c(1,0))
    
    dxs <- sp::spDistsN1(par.df[,poscoords], par.df[i,poscoords], longlat = longlat)
    if (!is.finite(dxs[i]))
        dxs[i] <- .Machine$double.xmax/2
  
    if(kernel=="gaussian")
        w.i<-weight.gaussian(dxs, bandwidth)
  
    if(kernel=="bisquare")
        w.i<-weight.bisquare(dxs, bandwidth)

    w.i[i]<-0
    w.i<-w.i*par.df[,posweights]
  
    if (any(w.i < 0 | is.na(w.i)))
        stop(paste("Invalid weights for i:", i))
    lm.i <- try(lm.wfit(y = par.df[,posy], x = par.df[,posx], w = w.i))
    
    if (!inherits(lm.i, "try-error")) {
        b <- coefficients(lm.i)
        return(par.df[i,posweights] * par.df[i,posy] - (t(b) %*% (par.df[i,posweights] * 
                                                            par.df[i, posx])))
  }
}
