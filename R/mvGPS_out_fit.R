#' Fit linear outcome regression using potential weights
#'
#' @inheritParams mvGPSsim_bal
#' @param Y numeric vector specifying the outcome of interest
#' @param W list of weights to use in the outcome 
#' @inheritParams hull_sample
#' @param poly_degree integer scalar specifying the polynomial degree which will 
#' be applied to each exposure and included in the outcome modeling. Default is NULL
#' signifying that only first order terms are included
#' @param D_interact logical indicator denoting whether to include all pairwise
#' interactions of univariate exposures when modeling
#'
#' @details Availables models include "mvGPS", "entropy", "CBPS", "PS"
#' 
#' @importFrom stats poly as.formula
#'
#' @export
mvGPS_out_fit <- function(Y, W, D, C, num_grid_pts=500, poly_degree=NULL, D_interact=FALSE){
    perf_metrics <- match.arg(perf_metrics, c("MSE", "bias"), several.ok=TRUE)
    if(!is.list(W)) stop("`W` must be a list with each element corresponding to weights to be used in outcome regression", call.=FALSE)
    if(nrow(D)!=nrow(C)) stop("`D` and `C` should have the same number of rows which corresponds to the number of units", call.=FALSE)
    n <- nrow(D)
    p <- ncol(D)
    #this is incase the object W is a data.frame or tibble and we want to convert it directly into a list
    W <- as.list(W)
    W_length <- unlist(lapply(W, length))
    if(!all(W_length==n)) stop("All weights in `W` must be same length as number of units in exposure and confounders. Check length of elements in `W`.", call.=FALSE)
    if(length(Y)!=n) stop("Vector `Y` must have same number of units as `D` and `C`", call.=FALSE)
    
    #calculating the convex hull data points
    hull_results <- hull_sample(D, num_grid_pts)
    hull_grid_pts <- hull_results$grid_pts
    colnames(hull_grid_pts) <- paste0("D", seq_len(p))
    
    #polynomials only
    if(!is.null(poly_degree) && !D_interact){
        #transforming the exposure
        D_poly <- lapply(seq_len(p), function(x) poly(D[, x], degree=poly_degree, raw=TRUE, simple=TRUE))
        D_poly <- do.call(cbind, D_poly)
        colnames(D_poly) <- paste0("D", rep(seq_len(p), each=poly_degree), paste0("_poly_", seq_len(poly_degree)))
        D <- D_poly
        
        #transforming the grid points
        hull_poly <- lapply(seq_len(p), function(x) poly(hull_grid_pts[, x], degree=poly_degree, raw=TRUE, simple=TRUE))
        hull_poly <- do.call(cbind, hull_poly)
        colnames(hull_poly) <- paste0("D", rep(seq_len(p), each=poly_degree), paste0("_poly_", seq_len(poly_degree)))
        hull_grid_pts <- hull_poly
    }
    #interaction of first order terms only
    if(D_interact==TRUE && is.null(poly_degree)){
        #transforming exposure
        D_int <- model.matrix(as.formula(paste0("~(", paste(colnames(D), collapse="+"), ")^2-1")), data=data.frame(D))
        D <- D_int
        
        #transforming grid points
        hull_int <- model.matrix(as.formula(paste0("~(", paste(colnames(hull_grid_pts), collapse="+"), ")^2-1")), data=data.frame(hull_grid_pts))
        hull_grid_pts <- hull_int
    }
    #interactions with polynomials
    if(D_interact==TRUE && !is.null(poly_degree)){
        D_poly_int <- model.matrix(as.formula(paste0("~", paste(paste0("poly(",colnames(D), ", degree=", poly_degree, ", raw=TRUE, simple=TRUE)"), collapse="*"), "-1")), data=data.frame(D))
        D <- D_poly_int
        
        #transforming grid points
        hull_poly_int <- model.matrix(as.formula(paste0("~", paste(paste0("poly(",colnames(hull_grid_pts), ", degree=", poly_degree, ", raw=TRUE, simple=TRUE)"), collapse="*"), "-1")), data=data.frame(hull_grid_pts))
        hull_grid_pts <- hull_poly_int
    }
    
    #constructing the design matrix
    X <- cbind(intcpt=rep(1, n), C, D)
    
    #fitting an unadjusted model which accounts does not have any weights or control for confounders
    unadj_mod <- lm(Y~D)
    unadj_beta_hat <- coef(unadj_mod)[-1]
    unadj_y_hat <- as.numeric(hull_grid_pts %*% unadj_beta_hat)
    
    mod_fit <- lapply(W, function(w){
        w_mod <- lm(Y~D, weights=w)
        beta_hat <- coef(w_mod)[-1]
        y_hat <- as.numeric(hull_grid_pts %*% beta_hat)
        return(list(beta_hat=beta_hat, y_hat=y_hat))
    })
    
    mod_fit[["unadj"]] <- list(beta_hat=unadj_beta_hat, y_hat=unadj_y_hat)
    
    return(list(mod_fit=mod_fit, hull_grid_pts=hull_grid_pts))
}
