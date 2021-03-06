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
#' @param collinearity Logical. If TRUE, metrics of local collinearity (local correlation and local VIF are computed)
#' @return A fitted GWR
#' @import snowfall
#' @import reshape2
#' @importFrom snow setDefaultClusterOptions
#' @import sp
#' @export
gwr_par<-function(formula, data, coords, bandwidth, weights=NULL,
	kernel="gaussian", longlat=NULL, se.fit=FALSE, diagnostic=FALSE, 
	collinearity=FALSE, adapt=F, ncores=NULL)
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
		    if(!is.na(is.projected(data)) && !is.projected(data))
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
    
    if(collinearity && ncol(x)<=2){
        warning("Number of variables < 2. Collinearity will not be computed")
        collinearity<-FALSE
    }
    
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
        	                            longlat=longlat, adapt=adapt, se.fit=se.fit, 
        	                            diagnostic=diagnostic, collinearity=collinearity))

    # Running linear models sequentially if ncores>2
    if(!is.null(ncores) && ncores>1){
    	snowfall::sfInit(cpus=ncores, parallel=TRUE)
    	snowfall::sfExport(list=c("x", "y", "coords", "bandwidth",
    		"weights", "kernel", "longlat", "adapt", "se.fit", "diagnostic"))
    	snowfall::sfLibrary(sp)
    	param.local.lm<-sfLapply(seq(n.sample), 
    		function(cell) gwr.internal(x=x, y=y, cell=cell, coords=coords, 
    		                            bandwidth=bandwidth, weights=weights,kernel=kernel,
    		                            longlat=longlat, adapt=adapt, se.fit=se.fit, 
    		                            diagnostic=diagnostic, collinearity=collinearity))
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
    tss<-c(cov.wt(matrix(y, ncol = 1), wt = weights, method = "ML")$cov * n.sample)
    #delta1<-sum(diag(crossprod2))
    #sigma2<-rss/delta1
    #odelta2<-sum(diag(crossprod2)^2)
    #delta2<-sum(diag(crossprod2 %*% crossprod2))
    sigma2.b<-rss/n.sample
    AICc<-2 * n.sample * log(sqrt(sigma2.b)) + n.sample * 
    	log(2 * pi) + (n.sample * ((n.sample + diaghatmat)/(n.sample - 2 - diaghatmat)))
    AIC<-2 * n.sample * log(sqrt(sigma2.b)) + n.sample * log(2 * pi) + n.sample + diaghatmat

    diagnostics<-list(AIC=AIC, AICc=AICc, RSS=rss, TSS=tss, EDF=effective.df)
    }

    else
        diagnostics<-NULL
        
    local.R2<-sapply(seq(n.sample), function(cell) local.R2(y, cell, coords,
                                                            df[,"yhat"], longlat, adapt,
                                                            weights, kernel, bandwidth))
    df<-cbind(df, local.R2)
    
    # turning the df into a sdf projecting the final dataframe
    if(!is.null(projection))
        sdf<-SpatialPointsDataFrame(coords=coords, data=as.data.frame(df),
                                    proj4string=CRS(projection))
    
    else
        sdf<-SpatialPointsDataFrame(coords=coords, data=as.data.frame(df))      

    # computing collinearity metrics
    if(collinearity){
        list.cor.df<-lapply(param.local.lm, function(j) j[["loc.cor.i"]])
        list.vif.df<-lapply(param.local.lm, function(j) j[["vif.i"]])
        
        names.corr<-combn(colnames(x)[-1], 2, function(x) paste0("corr_", paste(x, collapse=".")))
        names.vif<-paste0("vif_", colnames(x)[-1])
        
        collin.df<-cbind(do.call(rbind, list.cor.df), do.call(rbind, list.vif.df))
        colnames(collin.df)<-c(names.corr, names.vif)
        
        rm(list=c("list.cor.df", "list.vif.df"))
        
        if(!is.null(projection))
            collin.sdf<-SpatialPointsDataFrame(coords=coords, data=as.data.frame(collin.df),
                                        proj4string=CRS(projection))
        else
            collin.sdf<-SpatialPointsDataFrame(coords=coords, data=as.data.frame(collin.df))
    }
    
    else
        collin.sdf<-NA
    
    results.list<-list(sdf=sdf, call=gwr.call, global.lm=lm.global, bandwidth=bandwidth,
                           kernel=kernel, diagnostic.metrics=diagnostics, collinearity=collin.sdf)
        
    class(results.list)<-"pargwr"
    invisible(results.list)
}

#'@export
print.pargwr<-function(x,...)
{
    if(class(x)!="pargwr")
        stop("Object class is not pargwr")
    
    cat("Call:\n")
    print(x$call)
    cat("Kernel:", x$kernel, "\n")
    cat("Bandwidth: ", x$bandwidth, "\n")
    
    coef_names<-c("Intercept", all.vars(formula(x$call))[-1])
    df_glm<-attr(logLik(x$global.lm), "df")
    n.sample<-attr(logLik(x$global.lm), "nobs")
    
    glob_param<-c(summary(x$global.lm)$coefficients[,"Estimate"],
                  summary(x$global.lm)$coefficients[,"Std. Error"],
                  as.integer(df_glm-1), as.integer(n.sample-df_glm), AIC(x$global.lm),
                  AIC(x$global.lm)+(2*df_glm^2+2*df_glm)/(n.sample-df_glm-1),
                  deviance(x$global.lm), summary(x$global.lm)$adj.r.squared)
    
    gwr_coeffs<-as(x$sdf, "data.frame")[, (1 + (1:(df_glm-1))), drop = FALSE]
    if(any(is.na(gwr_coeffs)))
        warning("NAs dropped in GWR local coefficients!")
    
    if(is.null(x$diagnostic.metrics))
        gwr_param<-c(apply(gwr_coeffs, 2, mean, na.rm=T), 
                     apply(gwr_coeffs, 2, sd, na.rm=T),
                     NA, NA, NA, NA, NA, NA)
    
    if(!is.null(x$diagnostic.metrics))
        gwr_param<-c(apply(gwr_coeffs, 2, mean, na.rm=T), 
                     apply(gwr_coeffs, 2, sd, na.rm=T),
                     x$diagnostic.metrics$EDF, n.sample-x$diagnostic.metrics$EDF,
                     x$diagnostic.metrics$AIC, x$diagnostic.metrics$AICc,
                     x$diagnostic.metrics$RSS, 
                     1-x$diagnostic.metrics$RSS/x$diagnostic.metrics$TSS)
    
    
    printDF<-rbind(glob_param, gwr_param)
    colnames(printDF)<-c(coef_names, paste0("sd_", coef_names), 
                         "Effective number of parameters",
                         "Effective DF", "AIC", "AICc",
                         "RSS", "(Quasi) R2")
    rownames(printDF)<-c("Global model", "GWR")
    
    printCoefmat(printDF)
    invisible(x)
}