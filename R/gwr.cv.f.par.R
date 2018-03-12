#### For any given bandwidth, gwr.cv.f.par compute the cv score (AIC score will come later as well as not fixed bandwidth)
# Done by running the n GLMs accross multiple cores/nodes.
# SHOULD NOT and cannot be run oustide from a call by gwr.sel.par()!

gwr.cv.f.par<-function (bandwidth, y, x, coords, kernel, verbose = TRUE, longlat = FALSE,
                        RMSE = FALSE, weights, show.error.messages = TRUE, ncores, cluster=cl1)
{
    n <- NROW(x)
    cv <- numeric(n)
    options(show.error.messages = show.error.messages)
    
    
    
    if(!is.null(ncores) && ncores>1){
        df<-cbind(y, x, weights, coords)
        desc.df<-bigmemory::describe(bigmemory::as.big.matrix(df))
        snow::clusterExport(cl=cluster, list=c("bandwidth", "kernel", "longlat", 
                                    "desc.df", "cv.compz", "weight.gaussian"),
                      envir=environment())
        snow::clusterEvalQ(cluster, library(bigmemory))
        
        cv<-snow::parSapply(cl=cluster, seq(n), function(m) cv.compz(m, descr=desc.df, longlat=longlat, 
                                                       kernel=kernel, bandwidth=bandwidth))
    }

    if(is.null(ncores) || ncores==1){
        desc.df<-cbind(y, x, weights, coords)
        cv<-sapply(seq(n), function(m) cv.compz(m, descr=desc.df, longlat=longlat, 
                                                kernel=kernel, bandwidth=bandwidth))
    }

    score <- sum(t(cv) %*% cv) #MSE
    if (RMSE)
        score <- sqrt(score/n) #RMSE
    if (!show.error.messages)
        options(show.error.messages = TRUE)
    if (verbose)
        cat("Bandwidth:", bandwidth, "CV score:", score, "\n")
    score
}