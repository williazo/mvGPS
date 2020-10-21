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
#' @param grid_type character value indicating the type of grid to sample from
#' the convex hull from \code{\link[sp]{spsample}}
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
#' exposure of interest. We can define the convex hull of these observations as
#' \strong{H}. There are two distinct processes for defining \strong{H} depending
#' on whether \eqn{m=2} or \eqn{m>2}.
#' 
#' If \eqn{m=2}, we use the \code{\link[grDevices]{chull}} function to define the 
#' vertices of the convex hull. The algorithm implemented is described in \insertCite{eddy1977new;textual}{mvGPS}.
#' 
#' If \eqn{m>2}, we use the \code{\link[geometry]{convhulln}} function. This algorithm 
#' for obtaining the convex hull in m-dimensional space uses Qhull described in
#' \insertCite{barber1996quickhull;textual}{mvGPS}. Currently this function returns
#' only the vertex set \code{hpts_vs} without the grid sample points. There are
#' options to visualize the convex hull when \eqn{m=3} using triangular facets,
#' but there are no implementable solutions to sample along convex hulls in higher
#' dimensions.
#' 
#' To restrict the convex hull to higher density regions of the exposure we can 
#' apply trimming. To apply trimming set \code{trim_hull=TRUE} and specify 
#' \code{trim_quantile=q} where \code{q} must be in \[0.5, 1\]. Along each
#' exposure dimension we then calculate the upper and lower bounds using the 
#' \code{\link[stats]{quantile}} function, i.e., \code{quantile(q)} and 
#' \code{quantile(1-q)}. Any observations that have a value above or below these
#' sample quantiles is excluded. The remaining observations that fall completely 
#' within the sample quantiles across all dimensions are used to estimate the 
#' convex hull. We return \code{X} that represents the observations used. 
#' If \code{trim_hull=FALSE}, then \code{X} is unchanged. However, if trimming
#' is applied then \code{X} contains only the remaining observations after trimming.
#' 
#' @examples 
#' #generating exposure with m=3
#' D <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(100))), nrow=100)
#' 
#' #first using only the first two variables we can return hpts_vs and grid_pts
#' D_hull <- hull_sample(D[, 1:2])
#' 
#' #when m>2 we only return hpts_vs and grid_pts is NULL
#' D_hull_large <- hull_sample(D)
#' is.null(D_hull_large$grid_pts)
#' 
#' #we can also apply trimming to the convex hull and return this reduced matrix
#' D_hull_trim <- hull_sample(D[, 1:2], trim_hull=TRUE, trim_quantile=0.95)
#' dim(D_hull$X)
#' dim(D_hull_trim$X)
#' 
#' #alternatively, we can also define the number of points to sample from for grid_pts
#' small_grid <- hull_sample(D[, 1:2], num_grid_pts=100)
#' length(D_hull$grid_pts)
#' length(small_grid$grid_pts)
#' 
#' @references
#'     \insertAllCited{}
#' 
#' @return 
#' \itemize{
#' \item \code{hpts_vs}: vertices of the convex hull in m-dimensional space
#' \item \code{grid_pts}: values of grid points sampled over the corresponding convex hull
#' \item \code{X}: data used to generate convex hull which may be trimmed
#' }
#' 
#' @export
hull_sample <- function(X,
                        num_grid_pts=500,
                        grid_type="regular",
                        trim_hull=FALSE,
                        trim_quantile=NULL) {
    X_rslt <- X_check(X)
    assign("X", X_rslt$X)
    assign("m", X_rslt$m)
    grid_type <- match.arg(grid_type, choices = c("regular", "random", "hexagonal"))
    if(trim_hull==TRUE){
        if(is.null(trim_quantile)) stop("trim_hull set to TRUE but trim_quantile not specified.", call.=FALSE)
        if(trim_quantile<0.5 | trim_quantile>1) stop("trim_quantile must be between [0.5, 1]", call.=FALSE)
        trim_upper <- apply(X, 2, quantile, trim_quantile)
        trim_lower <- apply(X, 2, quantile, 1 - trim_quantile)
        X_trim <- sapply(seq_len(m), function(x){
            ifelse(X[, x] > trim_upper[x], NA, 
                   ifelse(X[, x] < trim_lower[x], NA, X[, x]))
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
        sp_grid_pts <- spsample(sps, n=num_grid_pts, type=grid_type)
        grid_pts <- coordinates(sp_grid_pts)
        colnames(grid_pts) <- colnames(X)
    } else {
        hpts <- geometry::convhulln(X)
        #returns an n times m matrix where each row contains indices used to form facets
        hpts_ind <- unique(c(hpts))
        #selecting all of the unique row indices as this will form the vertex set
        hpts_vs <- X[hpts_ind, ]
        #subsetting the original data to form the vertex set
        #NOTE: this vertex set is not ordered like the case of m=2
        grid_pts <- NULL
        #need to investigate ways to sample along m-dimensional object formed using vertex set
    }
    return(list(hpts_vs=hpts_vs, grid_pts=grid_pts, X=X))
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