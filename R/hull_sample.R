#' Sample Points Along a Convex Hull
#' 
#' @param X numeric matrix of n by p dimensions. Each row corresponds to a point in p-dimensional space.
#' @param num_grid_pts integer scalar denoting the number of parameters to 
#' search for over the convex hull. Default is 500.
#' @param trim_hull logical indicator of whether to restrict convex hull. Default is FALSE
#' @param trim_quantile numeric scalar between \[0.5, 1\] representing the 
#' quantile value to trim the convex hull. Only used if trim_hull is set to TRUE.
#' 
#' @import geometry
#' @import sp
#' 
#' @importFrom grDevices chull
#' @importFrom stats na.omit
#' 
#' @export
hull_sample <- function(X,  num_grid_pts=500, trim_hull=FALSE, trim_quantile=NULL){
    requireNamespace("geometry")
    requireNamespace("sp")
    if(!is.matrix(X)) stop("`X` must be a numeric matrix", call.=FALSE)
    p <- ncol(X)
    if(trim_hull==TRUE){
        if(is.null(trim_quantile)) stop("trim_hull set to TRUE but trim_quantile not specified.", call.=FALSE)
        if(trim_quantile<0.5 | trim_quantile>1) stop("trim_quantile must be between [0.5, 1]", call.=FALSE)
        trim_upper <- apply(X, 2, quantile, trim_quantile)
        trim_lower <- apply(X, 2, quantile, 1 - trim_quantile)
        X <- sapply(seq_len(p), function(x){
            ifelse(X[, x]>trim_upper[x], NA, 
                   ifelse( X[, x]<trim_lower[x], NA, X[, x]))
        })
        X <- na.omit(X)
    }
    
    if(p==2){
        #special base operations to use if operating in two dimensional space
        hpts <-chull(X) #coordinates for the convex hull
        hpts <- c(hpts, hpts[1]) #closing the set
        hpts_vs <- as.matrix(X[hpts, ]) #vertices of the convex hull
        p <- Polygon(hpts_vs)
        # wrap Polygon into Polygons object
        ps <- Polygons(list(p), 1)
        
        # wrap Polygons object into SpatialPolygons object
        sps <- SpatialPolygons(list(ps))
        
        #sampling regular points along grid of polygon
        sp_grid_pts <- spsample(sps, n=num_grid_pts, type="regular")
        grid_pts <- coordinates(sp_grid_pts)
        
    } else if(p>2){
        hpts <- geometry::convhulln(X)
        ## need to come back and work on this process when dimension is greater than 2
    }
    return(list(hpts_vs=hpts_vs, grid_pts=grid_pts))
}