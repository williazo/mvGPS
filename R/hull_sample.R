#' Sample Points Along a Convex Hull
#' 
#' To define a proper estimable region with multivariate exposure we construct 
#' a convex hull of the data in order to maintain the positivity identifying assumption.
#' We also provide options to create trimmed versions of the convex hull to further
#' restrict to high density regions in multidimensional space.
#' 
#' @param X numeric matrix of n by m dimensions. Each row corresponds to a point in m-dimensional space.
#' @param num_grid_pts integer scalar denoting the number of parameters to 
#' search for over the convex hull. Default is 500.
#' @param trim_hull logical indicator of whether to restrict convex hull. Default is FALSE
#' @param trim_quantile numeric scalar between \[0.5, 1\] representing the 
#' quantile value to trim the convex hull. Only used if trim_hull is set to TRUE.
#' 
#' @import geometry
#' @import sp
#' @importFrom grDevices chull
#' @importFrom stats na.omit
#' 
#' @details 
#' Assume that \eqn{X} is an \eqn{n\times m} matrix representing the multivariate
#' exposure of interest. We can 
#' 
#' @export
hull_sample <- function(X,  num_grid_pts=500, trim_hull=FALSE, trim_quantile=NULL){
    X_rslt <- X_check(X)
    assign("X", X_rslt$X)
    assign("m", X_rslt$m)
    if(trim_hull==TRUE){
        if(is.null(trim_quantile)) stop("trim_hull set to TRUE but trim_quantile not specified.", call.=FALSE)
        if(trim_quantile<0.5 | trim_quantile>1) stop("trim_quantile must be between [0.5, 1]", call.=FALSE)
        trim_upper <- apply(X, 2, quantile, trim_quantile)
        trim_lower <- apply(X, 2, quantile, 1 - trim_quantile)
        X_trim <- sapply(seq_len(m), function(x){
            ifelse(X[, x]>trim_upper[x], NA, 
                   ifelse( X[, x]<trim_lower[x], NA, X[, x]))
        })
        colnames(X_trim) <- colnames(X)
        X <- na.omit(X_trim)
    }
    
    if(m==2){
        #special base operations to use if operating in two dimensional space
        hpts <-chull(X) #coordinates for the convex hull
        hpts <- c(hpts, hpts[1]) #closing the set
        hpts_vs <- as.matrix(X[hpts, ]) #vertices of the convex hull
        m <- Polygon(hpts_vs)
        # wrap Polygon into Polygons object
        ps <- Polygons(list(m), 1)
        
        # wrap Polygons object into SpatialPolygons object
        sps <- SpatialPolygons(list(ps))
        
        #sampling regular points along grid of polygon
        sp_grid_pts <- spsample(sps, n=num_grid_pts, type="regular")
        grid_pts <- coordinates(sp_grid_pts)
        
    } else {
        hpts <- geometry::convhulln(X)
        ## need to come back and work on this process when dimension is greater than 2
    }
    return(list(hpts_vs=hpts_vs, grid_pts=grid_pts))
}

#' Checking that the exposure matrix is properly specified
#' 
#' @inheritParams hull_sample
#' 
#' @keywords internal
X_check <- function(X){
    X <- as.matrix(X)
    m <- ncol(X)
    
    if(!any(apply(X, 2, is.numeric))) stop("X must be numeric", call.=FALSE)
    if(m<2) stop("Exposure is not multivariate", call.=FALSE)
    
    return(list(X=X, m=m))
}