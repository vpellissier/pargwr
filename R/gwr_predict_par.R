#' GWR prediction
#' 
#' This function predicts values of Y in location were it was not observed
#'   
#' @param fittedGWR A GWR object fitted with gwr_par()
#' @param newdata A new dataset used to predict the value (either data.frame of Spatial*Dataframe object)
#' @param newcoords A two-columns matrix with the coordinates of the new dataset if newdata is a dataframe Ignored otherwise
#' @param longlat TRUE if coordinates are longitude-latitude in decimal degrees, in which case, distances are measured in kilometers
#' @param bandwith Bandwidth. Can be computed with gwr.sel.par()
#' @param se.fit Logical. TRUE if standard errors of the fit should be assessed.
#' @param ncores Number of cores in which the computation should be spanned
#' @return A SDF containing predicted values
#' @import snowfall
#' @import sp
#' @export
gwr_predict_par<-function(fittedGWR=NULL, newdata, newcoords, 
	longlat=NULL, se.fit=FALSE, ncores=NULL)
{
	adapt<-FALSE # to be implemented later
	if(is.null(fittedGWR))
		stop("Could not find a fittedGWR object")

	if(!is(fittedGWR, "pargwr"))
		stop("fittedGWR is not an object of the class 'pargwr'")

	# Gathering data from the fittedGWR object
	x<-model.matrix(fittedGWR$global.lm)
	x.vars<-all.vars(formula(fittedGWR$call))[-1]
	y<-model.response(model.frame(fittedGWR$global.lm))
	weights<-model.extract(model.frame(fittedGWR$global.lm), "weights")

	coords<-coordinates(fittedGWR$sdf)
	bandwidth<-fittedGWR$bandwidth
	kernel<-fittedGWR$kernel
	yhat<-fittedGWR$sdf@data[,"yhat"]

	# Conducting check for the new dataset
    newproj<-NULL

	if(is(newdata, "Spatial")){
		if(!missing(newcoords))
			warning("Coordinates of the fitpoints are taken directly from data")
		newcoords<-coordinates(newdata)
		newproj<-proj4string(newdata)

		if(is.null(longlat) || !is.logical(longlat)){
			if(!is.na(is.projected(newdata) && !is.projected(newdata)))
				longlat<-TRUE
			else
				longlat<-FALSE			
		}

		newdata<-newdata@data
	}

	if (missing(newcoords))
		stop("Please provide a the coordinates of the fitpoints (newcoords)")

	if(is.null(longlat) || !is.logical(longlat))
		longlat<-FALSE

	newdata<-try(newdata[,x.vars, drop=FALSE])
	if(class(newdata)=="try-error")
		stop("Names of the new dataset do not match the names in the fittedGWR object")

	n.new<-nrow(newdata)
	newdata<-cbind(rep(1, n.new), newdata)
	colnames(newdata)[1]<-"(Intercept)"

    # Running linear models sequentially if ncores==NULL
    if(is.null(ncores))
        param.local.lm<-lapply(seq(n.new), 
        	function(cell) gwr.pred.internal(x, y, coords, cell, newcoords,
    			newdata, bandwidth, weights, yhat,
    			kernel, longlat, adapt, se.fit))

    # Running linear models sequentially if ncores>2
    if(!is.null(ncores) && ncores>1){
    	snowfall::sfInit(cpus=ncores, parallel=TRUE)
    	snowfall::sfExport(list=c("x", "y", "coords", "cell", "newcoords",
    			"newdata", "bandwidth", "weights", "yhat",
    			"kernel", "longlat", "adapt", "se.fit"))
    	snowfall::sfLibrary(sp)
    	param.local.lm<-sfLapply(seq(n.new), 
    		function(cell) gwr.pred.internal(x, y, coords, cell, newcoords,
    			newdata, bandwidth, weights, yhat,
    			kernel, longlat, adapt, se.fit))
    	snowfall::sfStop()
    }

    list.df<-lapply(param.local.lm, function(j) j[["df.i"]])
    df<-do.call(rbind, list.df)

    if(!is.null(newproj))
        sdf<-SpatialPointsDataFrame(coords=newcoords, data=as.data.frame(df),
    		proj4string=CRS(newproj))
    else
        sdf<-SpatialPointsDataFrame(coords=newcoords, data=as.data.frame(df))      

   	results.list<-list(sdf=sdf, bandwidth=bandwidth,
    	kernel=kernel)
    class(results.list)<-"pargwr"
    invisible(results.list)
}