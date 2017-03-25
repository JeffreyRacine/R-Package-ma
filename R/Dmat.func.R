Dmat.func <-
function(model,method=c("mma","jma")) {
    method <- match.arg(method)
    ## Residuals (Mallows model averaging)
    if(method=="mma") return(residuals(model))
    ## Jackknife fitted values (jackknife model averaging)
    htt <- hatvalues(model)
    if(method=="jma") return(fitted(model) - htt*residuals(model)/(1-htt))    
}
