#' GWR computation
#' 
#' This function computes a GWR (Geographically Weighted Regression). It spans the computation across several nodes in a cluster.
#'  
#' @param formula Formula of the GWR
#' @param data Dataset (either data.frame of Spatial*Dataframe object)
#' @param coords A two-columns matrix with the coordinates as X-Y if data is a data.frame
#' @param adapt Logical. TRUE if Adaptative bandwith, FALSE if fixed
#' @param kernel Character string. Weight kernel. Either "gaussian", "bisquare" or "boxcar"
#' @param longlat TRUE if coordinates are longitude-latitude in decimal degrees, in which case, distances are measured in kilometers
#' @param bandwith Bandwidth. Can be computed with gwr.sel.par()
#' @param se.fit Logical. TRUE if standard errors of the fit should be assessed.
#' @param ncores Number of cores in which the computation should be spanned
#' @param diagnostic Logical. If TRUE, the following metrics are displayed: AIC, AICc, RSS, Effective numbers of parameters and of degrees of freedom.
#' @return A fitted GWR
#' @import snowfall
#' @importFrom snow setDefaultClusterOptions
#' @import sp
#' @export
gwr_par<-function(formula, data, coords, bandwidth, weights=NULL,
	kernel="gaussian", longlat=NULL, se.fit=FALSE, diagnostic=FALSE, adapt=F, ncores=NULL)
{
	if(missing(formula)) 
		stop("Formula is missing")
	
	if (missing(bandwidth))
		stop("Please provide a bandwidth. Can be computed using gwr.sel.par()")

    gwr.call<-match.call()
    projection<-NULL

	if(is(data, "Spatial")){
		if(!missing(coords))
			warning("Coordinates are taken directly from data")
		coords<-coordinates(data)
		projection<-proj4string(data)

		if(is.null(longlat) || !is.logical(longlat)){
			if(!is.na(is.projected(data) && !is.projected(data)))
				longlat<-TRUE
			else
				longlat<-FALSE			
		}

		data<-data@data
	}

	if (missing(coords))
		stop("Please provide a coordinates matrix")

	if(is.null(longlat) || !is.logical(longlat))
		longlat<-FALSE

	y.var<-all.vars(formula)[1]
	x.vars<-all.vars(formula)[-1]
	n.vars<-length(x.vars)

	y<-as.vector(data[,y.var])
	n.sample<-length(y)

	x<-data[,x.vars, drop=FALSE]
    x<-cbind(rep(1, n.sample), x)
    colnames(x)[1]<-"(Intercept)"
    x<-as.matrix(x)
    
	if(!is.null(weights) && !is.numeric(weights))
		stop("Weights should be a numeric vector")

	if(is.null(weights))
		weights<-rep(1, n.sample)

	lm.global<-lm(y~x[,-1], weights=weights)

	# Settings to run local linear models
	coeffs<-matrix(nrow=length(y), ncol=ncol(x))
    colnames(coeffs)<-colnames(x)

    # The hatmatrix needs to be computed if diagnostic=TRUE (AICc, AIC, sigma...)

    # Running linear models sequentially if ncores==NULL
    if(is.null(ncores))
        param.local.lm<-lapply(seq(n.sample), 
        	function(cell) gwr.internal(x=x, y=y, cell=cell, coords=coords, 
    			bandwidth=bandwidth, weights=weights,kernel=kernel, 
    			longlat=longlat, adapt=adapt, se.fit=se.fit, diagnostic=diagnostic))

    # Running linear models sequentially if ncores>2
    if(!is.null(ncores) && ncores>1){
    	snowfall::sfInit(cpus=ncores, parallel=TRUE)
    	snowfall::sfExport(list=c("x", "y", "coords", "bandwidth",
    		"weights", "kernel", "longlat", "adapt", "se.fit", "diagnostic"))
    	snowfall::sfLibrary(sp)
    	param.local.lm<-sfLapply(seq(n.sample), 
    		function(cell) gwr.internal(x=x, y=y, cell=cell, coords=coords, 
    			bandwidth=bandwidth, weights=weights,kernel=kernel, 
    			longlat=longlat, adapt=adapt, se.fit=se.fit, diagnostic=diagnostic))
    	snowfall::sfStop()
    }

    list.df<-lapply(param.local.lm, function(j) j[["df.i"]])
    df<-do.call(rbind, list.df)

    if(diagnostic==TRUE){
        list.hatmat<-lapply(param.local.lm, function(j) j[["lhat.i"]])
        hatmat<-do.call(rbind, list.hatmat)
    
    # This bloc computes diagnostic metrics, based on the hatmatrix
    diaghatmat<-sum(diag(hatmat))
    crossprod1<-t(hatmat) %*% hatmat
    diagcrossprod<-sum(diag(crossprod1))
    effective.df<-n.sample - 2 * diaghatmat + diagcrossprod
    crossprod2<-t(diag(n.sample) - hatmat) %*% (diag(n.sample) - hatmat)
    rss<-c(t(y) %*% crossprod2 %*% y)
    #delta1<-sum(diag(crossprod2))
    #sigma2<-rss/delta1
    #odelta2<-sum(diag(crossprod2)^2)
    #delta2<-sum(diag(crossprod2 %*% crossprod2))
    sigma2.b<-rss/n.sample
    AICc<-2 * n.sample * log(sqrt(sigma2.b)) + n.sample * 
    	log(2 * pi) + (n.sample * ((n.sample + diaghatmat)/(n.sample - 2 - diaghatmat)))
    AIC<-2 * n.sample * log(sqrt(sigma2.b)) + n.sample * log(2 * 
        pi) + n.sample + diaghatmat

    diagnostics<-list(AIC=AIC, AICc=AICc, RSS=rss, EDF=effective.df)
    }

    else
        diagnostics<-NULL

    local.R2<-sapply(seq(n.sample), function(cell) local.R2(y, cell, coords, 
    													df[,"yhat"], longlat, adapt,
    													weights, kernel, bandwidth))
    df<-cbind(df, local.R2)
    
    if(!is.null(projection))
        sdf<-SpatialPointsDataFrame(coords=coords, data=as.data.frame(df),
    		proj4string=CRS(projection))
    else
        sdf<-SpatialPointsDataFrame(coords=coords, data=as.data.frame(df))      

   	results.list<-list(sdf=sdf, call=gwr.call, global.lm=lm.global, bandwidth=bandwidth,
    	kernel=kernel, diagnostic.metrics=diagnostics)
    class(results.list)<-"pargwr"
    invisible(results.list)
}