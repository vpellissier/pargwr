#' GWR computation
#' 
#' This function computes a GWR (Geographically Weighted Regression). It spans the computation across several nodes in a cluster.
#'  
#' @param formula Formula of the GWR
#' @param data Dataset (either data.frame of SpatialDataframe object)
#' @param coords A two-columns matrix with the coordinates as X-Y if data is a data.frame
#' @param adapt Logical. TRUE if Adaptative bandwith, FALSE if fixed
#' @param kernel Character string. Weight kernel. Either "gaussian", "bisquare" or "boxcar"
#' @param longlat TRUE if coordinates are longitude-latitude in decimal degrees, in which case, distances are measured in kilometers
#' @param bandwith Minimum distance between two separate runs of the optimize
#' @param se.fit Logical. TRUE if standard errors of the fit should be assesse. Not implemente yet
#' @param ncores Number of cores in which the computation should be spanned
#' @return A fitted GWR
#' @export
gwr_par<-function(formula, data, coords, bandwidth, weights=NULL,
	kernel="gaussian", longlat=NULL, se.fit=FALSE, adapt=F, ncores=NULL)
{
	if(missing(formula)) 
		stop("Formula is missing")
	
	if (missing(bandwidth))
		stop("Please provide a bandwidth. Can be computed using gwr.sel.par()")

	if(is(data, "Spatial")){
		if(!missing(coords))
			warning("Coordinates are taken directly from data")
		coords<-coordinates(data)
		projection<-proj4string(data)

		if(is.null(longlat) || !is.logical(longlat)){
			if(!is.projected(data))
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
	
	y<-as.vector(data[,y.var])
	n.sample<-length(y)

	x<-as.matrix(data[,x.vars])
	Intercept<-rep(1, n.sample)
	x<-cbind(Intercept, x)

	if(!is.null(weights) && !is.numeric(weights))
		stop("Weights should be a numeric vector")

	if(is.null(weights))
		weights<-rep(1, n.sample)

	lm.global<-lm.wfit(x, y, w=weights)

	# Settings to run local linear models
	coeffs<-matrix(nrow=length(y), ncol=ncol(x))
    colnames(coeffs)<-colnames(x)

    # Running linear models sequentially if ncores==NULL
    if(is.null(ncores))
        param.local.lm<-lapply(seq(n.sample), function(cell) gwr.internal(x, y, cell,
    														coords, bandwidth, weights, 
    														kernel, longlat, adapt, se.fit))
    # Running linear models sequentially if ncores>2
    if(!is.null(ncores) && ncores>1){
    	snowfall::sfInit(cpus=ncores, parallel=TRUE)
    	snowfall::sfExport(list=c("x", "y", "coords", "bandwidth",
    		"weights", "kernel", "longlat", "adapt", "se.fit"))
    	snowfall::sfLibrary(sp)
    	param.local.lm<-sfLapply(seq(n.sample), function(cell) gwr.internal(x, y, cell,
    														coords, bandwidth, weights, 
    														kernel, longlat, adapt, se.fit))
    	snowfall::sfStop()
    }

    df<-do.call(rbind, param.local.lm)

    local.R2<-sapply(seq(n.sample), function(cell) local_R2(y, cell, coords, 
    													df[,"yhat"], longlat, adapt,
    													weights, kernel, bandwidth))
    df<-cbind(df, local.R2)

    return(df)
}