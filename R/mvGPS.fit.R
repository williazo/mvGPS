#' Fit Weighted Linear Regression for Dose-Response
#' 
#' Fitting the dose-response model for causal inference model using weights
#' generated from various propensity score models. Can use either the 
#' \code{\link[mvGPS]{mvGPS}} function to generate weights or the 
#' \code{\link[mvGPS]{bal}} function to generate entire list of weights from
#' different propensity models.
#'
#' @param formula  object of class \code{formula} which symbolically describes 
#' the dose-response with the outcome on one side and a function of the exposures on
#' the right side.
#' @param W list of weights to use in the dose-response regression
#' @param data data.frame that contains the outcome and covariates for 
#' dose-response. If none is provided then assumes objects are in user's environment.
#'
#' @details 
#' Will describe details here
#' 
#' @importFrom stats poly as.formula
#'
#' @export
mvGPS.fit <- function(formula, W, data = sys.frame(sys.parent())){
    #checks for proper specification of formula
    if(class(formula) != "formula"){
        stop("formula improperly specified. See details", call. = FALSE)
    } else if (length(formula) != 3){
        stop("formula must be 2-sided", call. = FALSE)
    }
    #outcome vector
    Y <- model.frame(formula, data)[, 1]
    #design matrix
    X <- model.matrix(formula, data)
    n <- length(Y)
    
    W <- as.list(W)
    W_length <- unlist(lapply(W, length))
    if(!all(W_length==n)) stop("all weights in W must be same length as number of units in exposure and confounders", call.=FALSE)

    #fitting an unadjusted model which accounts does not have any weights or control for confounders
    unadj_mod <- lm(Y ~ X - 1)
    
    mod_fit <- lapply(W, function(wts) lm(Y ~ X - 1, weights=wts))
    # mod_fit <- lapply(W, function(wts) gls(model=formula, data=mf, weights=varFunc(~wts)))
    #adding the unadjusted model
    mod_fit[["unweighted"]] <- list(mod=unadj_mod)
    
    mod_fit
}
