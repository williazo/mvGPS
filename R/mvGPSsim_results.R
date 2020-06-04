#' Construct Summary Statistics for Simulation Models
#'
#' @param D numeric matrix of dimension \code{n} by \code{m} designating values of the exposure
#' @param C numeric matrix of dimension \code{n} by \code{k} designating values of the confounders for each exposure value
#' @param W list of weights to use in the outcome 
#' @param alpha numeric vector of length \code{k+2+1} used for constructing the mean of the outcome \code{Y}.
#' The first term represents the intercept, the next \code{k} terms denote the effect of the confounders,
#' and the final \code{2} terms represent the true treatment effect
#' where \code{k} represents the number of confounders
#' @param perf_metrics character vector of which performance metrics to return. 
#' Default is to return "MSE" , "hull MSE", and "bias". See details for other options.
#' @param sd_Y numeric scalar designating the standard deviation to assume when
#' generating the outcome. Error term is assummed to be normally distributed.
#' @inheritParams hull_sample
#' @param poly_degree integer scalar specifying the polynomial degree which will 
#' be applied to each exposure and included in the outcome modeling. Default is NULL
#' signifying that only first order terms are included
#' @param D_interact logical indicator denoting whether to include all pairwise
#' interactions of univariate exposures when modeling
#'
#' @details Availables models include "mvGPS", "entropy", "CBPS", "PS"
#' 
#' @importFrom stats poly as.formula predict
#'
#' @export
mvGPSsim_results <- function(W, D, C, alpha, sd_Y, num_grid_pts=500, trim_hull=FALSE, trim_quantile=NULL,
                             perf_metrics=c("MSE", "bias", "hull MSE"), poly_degree=NULL, D_interact=FALSE){
    perf_metrics <- match.arg(perf_metrics, c("MSE", "bias",  "hull MSE"), several.ok=TRUE)
    n <- nrow(D)
    m <- ncol(D)
    C_n <- nrow(C)
    k <- ncol(C)
    if(!is.list(W)) stop("`W` must be a list with each element corresponding to weights to be used in outcome regression", call.=FALSE)
    if(m<2) stop("'D' must be exposure matrix with number of columns, m, greater than or equal to 2", call.=FALSE)
    if(!all(C_n==n)) stop("Each matrix in `C` must have same number of observations `n` as `D`", call.=FALSE)
    
    #this is in case the object W is a data.frame or tibble and we want to convert it directly into a list
    W <- as.list(W)
    W_length <- unlist(lapply(W, length))
    if(!all(W_length==n)) stop("All weights in `W` must be same length as number of units in exposure and confounders. Check length of elements in `W`.", call.=FALSE)
    
    
    #calculating the convex hull data points
    hull_results <- hull_sample(D, num_grid_pts, trim_hull, trim_quantile)
    hull_grid_pts <- hull_results$grid_pts
    colnames(hull_grid_pts) <- paste0("D", seq_len(m))
    
    #polynomials only
    if(!is.null(poly_degree) && !D_interact){
        #transforming the exposure
        D_poly <- lapply(seq_len(m), function(x) poly(D[, x], degree=poly_degree, raw=TRUE, simple=TRUE))
        D_poly <- do.call(cbind, D_poly)
        colnames(D_poly) <- paste0("D", rep(seq_len(m), each=poly_degree), paste0("_poly_", seq_len(poly_degree)))
        D <- D_poly
        
        #transforming the grid points
        hull_poly <- lapply(seq_len(m), function(x) poly(hull_grid_pts[, x], degree=poly_degree, raw=TRUE, simple=TRUE))
        hull_poly <- do.call(cbind, hull_poly)
        colnames(hull_poly) <- paste0("D", rep(seq_len(m), each=poly_degree), paste0("_poly_", seq_len(poly_degree)))
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
    m <- ncol(D)
    if(length(alpha)!=(1+k+m)){
        stop(paste0("`alpha` is not correctly specified. Needs to be length ", 1+k+m, "\n
             Transformed `D` matrix has ", m, " columns"), call.=FALSE)
        }
    alpha_D <- alpha[(2+k):length(alpha)]
    #generating the outcome Y as a linear model of C and D given coefficients alpha
    X <- cbind(intcpt=rep(1, n), C, D)
    Y <- X%*%alpha + rnorm(n, sd=sd_Y)
    
    #generating the ground truth for the points along the sample grid
    true_Y <- as.numeric(hull_grid_pts%*%alpha_D)
    
    #fitting an unadjusted model which accounts does not have any weights or control for confounders
    unadj_mod <- lm(Y~D)
    unadj_beta_hat <- coef(unadj_mod)[-1]
    unadj_y_hat <- predict(unadj_mod)
    hull_unadj_y_hat <- as.numeric(hull_grid_pts %*% unadj_beta_hat)
    
    mod_fit <- lapply(W, function(w){
        w_mod <- lm(Y~D, weights=w)
        beta_hat <- coef(w_mod)[-1]
        y_hat <- predict(w_mod)
        hull_y_hat <- as.numeric(hull_grid_pts %*% beta_hat)
        return(list(beta_hat=beta_hat, y_hat=y_hat, hull_y_hat=hull_y_hat))
    })
    
    mod_fit[["unadj"]] <- list(beta_hat=unadj_beta_hat, 
                               y_hat=unadj_y_hat, 
                               hull_y_hat=hull_unadj_y_hat)
    
    if("bias"%in%perf_metrics){
        bias_stats <- lapply(mod_fit, function(m){
            as.numeric(alpha_D - m$beta_hat)
        })
        bias_stats <- do.call(rbind, bias_stats)
        colnames(bias_stats) <- paste0(colnames(D), "_bias")
    } else{
        bias_stats <- NULL
    }
    
    if("hull MSE"%in%perf_metrics){
        hull_mse_stats <- lapply(mod_fit, function(m){
            as.numeric(crossprod(true_Y-m$hull_y_hat)/length(m$hull_y_hat))
        })
        hull_mse_stats <- unlist(hull_mse_stats)
    } else {
        hull_mse_stats <- NULL
    }
    
    if("MSE"%in%perf_metrics){
        mse_stats <- lapply(mod_fit, function(m){
            as.numeric(crossprod(Y-m$y_hat)/length(m$y_hat))
        })
        mse_stats <- unlist(mse_stats)
    } else {
        mse_stats <- NULL
    }
    return(list(bias_stats=bias_stats, mse_stats=mse_stats, hull_mse_stats=hull_mse_stats))
}
