#' GWR bandwidth selection
#' 
#' This function computes the optimal bandwidth for a GWR (Geographically Weighted Regression). It spans the computation across several nodes in a cluster.
#'  
#' @param formula Formula of the GWR
#' @param data Dataset (either data.frame of SpatialDataframe object)
#' @param coords A two-columns matrix with the coordinates as X-Y if data is a data.frame
#' @param adaptative Logical. TRUE if adaptative bandwith (k nearest neighbors), FALSE if fixed
#' @param kernel Character string. Weight kernel. Either "gaussian" or "bisquare"
#' @param method Validation method. So far, only the cross-validation approach is implemented 
#' @param longlat TRUE if coordinates are longitude-latitude in decimal degrees, in which case, distances are measured in kilometers
#' @param interval_dist Minimum distance between two separate runs of the optimizer (meters if adaptative=F; number of neibhbors if adaptative=T).
#' @param min_dist Minimum bandwith (meters if adaptative=F; number of neibhbors if adaptative=T). 
#' @param max_dist Maximum bandwith (meters if adaptative=F; number of neibhbors if adaptative=T)
#' @return A bandwidth
#' @export
gwr_sel_par<-function (formula, data = list(), coords, adaptative = FALSE, kernel="gaussian", 
                       method = "cv", verbose = TRUE, longlat = NULL, RMSE = FALSE, 
                       weights=NULL, interval_dist = 100, show.error.messages = TRUE, 
                       ncores = NULL, min_dist=NULL, max_dist=NULL) 
{
    if(missing(formula)) 
        stop("Formula is missing")

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
    
    if(!is.null(weights) && !is.numeric(weights))
        stop("Weights should be a numeric vector")
    
    if(is.null(weights))
        weights<-rep(1, n.sample)
            
    if (!adaptative){
        bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
        difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
        if (any(!is.finite(difmin))) 
            difmin[which(!is.finite(difmin))] <- 0
        if(is.null(min_dist)) 
            min_dist <- difmin/1000
        if (is.null(max_dist)) 
            max_dist <- difmin
        }
    
    if(adaptative){
        if(is.null(max_dist))
            max_dist <- length(y)
        if(is.null(min_dist))
            min_dist<-1
        }  
        
    if(is.null(ncores) || ncores==1){
            opt <- optimize(gwr.cv.f.par, lower=min_dist,upper=max_dist, 
                            maximum = FALSE, y = y, x = x, coords = coords,
                            adaptative = adaptative, kernel = kernel, verbose = verbose, 
                            longlat = longlat, RMSE = RMSE, weights = weights, 
                            ncores=ncores, show.error.messages = show.error.messages, 
                            tol = interval_dist*3)
            }
    
        if(!is.null(ncores) && ncores>1){    
            snowfall::sfInit(parallel=TRUE, cpus=ncores)
            snowfall::sfExport(list=c("coords", "longlat", "x", "y", "weights", "kernel"))
            snowfall::sfLibrary(sp)
    
            opt <- optimize(gwr.cv.f.par, lower=min_dist,upper=max_dist, 
                            maximum = FALSE, y = y, x = x, coords = coords,
                            adaptative = adaptative, kernel = kernel, verbose = verbose, 
                            longlat = longlat, RMSE = RMSE, weights = weights, 
                            ncores=ncores, show.error.messages = show.error.messages, 
                            tol = interval_dist*3)
    
            snowfall::sfStop()
            }

        res<-opt$minimum
        res
}