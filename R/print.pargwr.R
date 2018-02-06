# Print function used to display basics informations when a pargwr object is called

function(x,...)
{
    if(class(x)!="pargwr")
        stop("Object class is not pargwr")
    
    cat("Call:\n")
    print(x$call)
    cat("Kernel:", x$kernel)
    cat("Bandwidth: ", x$bandwidth)
    
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
    
    if(is.null(gwr2$diagnostic.metrics))
        gwr_param<-c(apply(gwr_coeffs, 2, mean, na.rm=T), 
                     apply(gwr_coeffs, 2, sd, na.rm=T),
                     NA, NA, NA, NA, NA, NA)
    
    if(!is.null(gwr2$diagnostic.metrics))
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
    
    printCoefmat(print(DF))
    invisible(x)
}