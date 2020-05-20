#' Sample Points Along a Convex Hull
#' 
#' @param X numeric matrix of n by p dimensions. Each row corresponds to a point in p-dimensional space.
#' @param num_grid_pts integer scalar denoting the number of parameters to 
#' search for over the convex hull. Default is 500.
#' 
#' @import geometry
#' @import sp
#' 
#' @importFrom grDevices chull
#' 
#' @export
hull_sample <- function(X,  num_grid_pts=500){
    requireNamespace("geometry")
    requireNamespace("sp")
    if(!is.matrix(X)) stop("`X` must be a numeric matrix", call.=FALSE)
    p <- ncol(X)
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