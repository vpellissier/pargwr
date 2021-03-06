#### For any given bandwidth, gwr.cv.f.par compute the cv score (AIC score will come later as well as not fixed bandwidth)
# Done by running the n GLMs accross multiple cores/nodes.
# SHOULD NOT and cannot be run oustide from a call by gwr.sel.par()!

gwr.cv.f.par<-function (bandwidth, y, x, coords, kernel, verbose = TRUE, longlat = FALSE,
                        RMSE = FALSE, weights, show.error.messages = TRUE, ncores)
{
    n <- NROW(x)
    cv <- numeric(n)
    options(show.error.messages = show.error.messages)

    if(!is.null(ncores) && ncores>1){
        snowfall::sfExport(list=c("bandwidth"))
        cv<-snowfall::sfSapply(seq(n), function(m) cv.compz(m, x=x, y=y, coords=coords, 
                                                            longlat=longlat, kernel=kernel, 
                                                            bandwidth=bandwidth, weights=weights))
    }

    if(is.null(ncores) || ncores==1){
        cv<-sapply(seq(n), function(m) cv.compz(m, x=x, y=y, coords=coords, 
                                                longlat=longlat, kernel=kernel, 
                                                bandwidth=bandwidth, weights=weights))
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