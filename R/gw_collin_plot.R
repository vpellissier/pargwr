#' Mapping of multicollinearity
#' 
#' Creates a map of multicollinearity indices (local correlation and local VIF) computed when the parameter 'collinearity' of the gwr_par function is set to TRUE.
#' 
#' @param obj A pargwr object fitted with gwr_par
#' @param metric String. The multicollinearity metric to be plotted (either 'vif' or 'corr')
#' @param plot Logical. Plot a map if TRUE, return a ggplot2 object if FALSE
#' @import ggplot2
#' @export
gw_collin_plot <- function(obj, metric = "vif"){
    if(class(obj)!="pargwr")
        stop("Object class is not pargwr")
    
    if(is.null(dim(obj$collinearity)))
        stop("Object does not contain collinearity metrics. Run gwr_par with 'collinearity'=TRUE")
    
    cutoff<-ifelse(metric == "vif", 10, 0.8)
    
    coords<-sp::coordinates(obj$collinearity)
    df<-obj$collinearity@data
    df<-df[,grep(metric, colnames(df))]
    names(df)<-gsub(paste0(metric, "_"), "", names(df))
    
    df.sp<-cbind(df, coords)
    df.sp<-reshape2::melt(df.sp, measure.vars=colnames(df))
    
    ggplot(df.sp, aes(x, y))+
            geom_raster(aes(fill=value>=cutoff))+
            facet_grid(~variable)
}
    