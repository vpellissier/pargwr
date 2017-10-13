#' GWR bandwidth selection
#' 
#' This function computes the optimal bandwidth for a GWR (Geographically Weighted Regression). It spans the computation across several nodes in a cluster.
#'  
#' @param formula Formula of the GWR
#' @param data Dataset (either data.frame of SpatialDataframe object)
#' @param coords A two-columns matrix with the coordinates as X-Y if data is a data.frame
#' @param adapt Logical. TRUE if Adaptative bandwith, FALSE if fixed
#' @param kernel Character string. Weight kernel. Either "gaussian" or "bisquare"
#' @param method Validation method. So far, only the cross-validation approach is implemented 
#' @param longlat TRUE if coordinates are longitude-latitude in decimal degrees, in which case, distances are measured in kilometers
#' @param interval_dist Minimum distance between two separate runs of the optimize
#' @param min_dist Minimum bandwith
#' @param max_dist Maximum bandwith
#' @return A bandwidth
#' @export
gwr_sel_par<-function (formula, data = list(), coords, adapt = FALSE, kernel="gaussian", 
                       method = "cv", verbose = TRUE, longlat = NULL, RMSE = FALSE, 
                       weights, interval_dist = 100, show.error.messages = TRUE, 
                       ncores = NULL, min_dist=NULL, max_dist=NULL) 
{
    if (!is.logical(adapt)) 
        stop("adapt must be logical")
    if (is(data, "Spatial")) {
        if (!missing(coords)) 
            warning("data is Spatial* object, ignoring coords argument")
        coords <- coordinates(data)
        if (is.null(longlat) || !is.logical(longlat)) {
            if (!is.na(is.projected(data)) && !is.projected(data)) {
                longlat <- TRUE
            }
            else {
                longlat <- FALSE
            }
        }
        data <- as(data, "data.frame")
    }
    if (is.null(longlat) || !is.logical(longlat)) 
        longlat <- FALSE
    if (missing(coords)) 
        stop("Observation coordinates have to be given")
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    dp.n <- length(model.extract(mf, "response"))
    weights <- as.vector(model.extract(mf, "weights"))
    
    if (is.null(weights)) 
        weights <- rep(as.numeric(1), dp.n)
    
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)   
    
    bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
    difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
    if (any(!is.finite(difmin))) 
        difmin[which(!is.finite(difmin))] <- 0
    if(is.null(min_dist)) min_dist <- difmin/1000
    if (is.null(max_dist)) max_dist <- difmin

    if(is.null(ncores)){
            opt <- optimize(gwr.cv.f.par, lower=min_dist,upper=max_dist, 
                    maximum = FALSE, y = y, x = x, coords = coords, 
                    kernel = kernel, verbose = verbose, longlat = longlat, 
                    RMSE = RMSE, weights = weights, show.error.messages = show.error.messages, 
                    tol = interval_dist*3)
    }
    
    if(!is.null(ncores) && ncores>1){    
    snowfall::sfInit(parallel=TRUE, cpus=ncores)
    snowfall::sfExport(list=c("coords", "longlat", "x", "y", "weights", "kernel"))
    snowfall::sfLibrary(sp)
    
    opt <- optimize(gwr.cv.f.par, lower=min_dist,upper=max_dist, 
                    maximum = FALSE, y = y, x = x, coords = coords, 
                    kernel = kernel, verbose = verbose, longlat = longlat, 
                    RMSE = RMSE, weights = weights, show.error.messages = show.error.messages, 
                    tol = interval_dist*3)
    
    snowfall::sfStop()
    }

    res<-opt$minimum
    res
}