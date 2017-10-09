#' GWR bandwidth selection
#' 
#' This function computes the optimal bandwidth for a given GWR using a cross-validation (leave-one-out) approach. It spans the computation across several nodes in a cluster.
#'  
#' @param formula Formula of the GWR
#' @param data Dataset (either data.frame of SpatialDataframe object)
#' @param coords A two-columns matrix with the coordinates as X-Y if data is a data.frame
#' @param adapt Logical. TRUE if Adaptative bandwith, FALSE if fixed
#' @param gweight Character string. Weight kernel. Either "gaussian" or "bisquare"
#' @return A bandwidth
#' @export
gwr.sel.par<-function (formula, data = list(), coords, adapt = FALSE, kernel="gaussian", 
                       method = "cv", verbose = TRUE, longlat = NULL, RMSE = FALSE, 
                       weights, tol = 300, show.error.messages = TRUE, 
                       ncores = 1, beta1=NULL, beta2=NULL) 
{
    require(snowfall, quietly=TRUE)
    require(spgwr, quietly=TRUE)
    
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
    if(is.null(beta1)) beta1 <- difmin/500  ## difmin = 537447 thus beta 1 = ca. 1km (1074.5m)
    if (is.null(beta2)) beta2 <- difmin / 50 ## == ca 10km (10748.9)
    
    
    
    sfInit(parallel=TRUE, cpus=ncores)
    sfExport(list=c("coords", "longlat", "x", "y", "weights", "gweight"))
    sfLibrary(sp)
    
    opt <- optimize(gwr.cv.f.par, lower=beta1,upper=beta2, 
                    maximum = FALSE, y = y, x = x, coords = coords, 
                    gweight = gweight, verbose = verbose, longlat = longlat, 
                    RMSE = RMSE, weights = weights, show.error.messages = show.error.messages, 
                    tol = tol)
    
    sfStop()
    res<-opt$minimum
    res
}