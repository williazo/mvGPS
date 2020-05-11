#' Construct Summary Statistics for Simulation Models
#'
#' @param model_list character string identifying which methods to use when constructing weights. See details for a list of appropriate models
#' @param all_uni logical indicator. If TRUE then all univariate models specified
#' in model_list will be estimated for each exposure. If FALSE will only estimate weights
#' for the first exposure
#' @param D numeric matrix of dimension \code{n} by \code{2} designating values of the exposure
#' @param C numeric matrix of dimension \code{n} by \code{k} designating values of the confounders
#' @param alpha numeric vector of length \code{k+2+1} used for constructing the mean of the outcome \code{Y}.
#' The first term represents the intercept, the next \code{k} terms denote the effect of the confounders,
#' and the final \code{2} terms represent the true treatment effect
#' where \code{k} represents the number of confounders
#'
#' @details Availables models include "mvGPS", "entropy", "CBPS", "PS"
#'
#' @export
mvGPSsim_results <- function(model_list, all_uni=TRUE, D, C, alpha){
    model_list <- match.arg(model_list, c("mvGPS", "entropy", "CBPS", "PS"), several.ok=TRUE)
}
