#' Multivariate Generalized Propensity Score
#'
#' Construct weights using propensity scores for multivariate continuous exposure
#' by assuming joint normal marginal and conditional densities. See details for
#' description of weights and additional methodology details.
#'
#' @param D numeric matrix of dimension \code{n} by \code{m} designating values of the exposures
#' @param C either a list of numeric matrices of length \code{m} of dimension 
#' \code{n} by \code{k} designating values of the confounders for each exposure 
#' value or if \arg{common} is TRUE a single matrix of of dimension \code{n} by
#' \code{k} that represents common confounders for all exposures.
#' @param trim_w logical indicator for whether to trim weights. default is FALSE
#' @param trim_quantile numeric scalar used to specify the upper quantile to 
#' trim weights if applicable. default is 0.99
#'
#' @importFrom stats model.frame model.matrix lm coef dnorm quantile
#' 
#' @return list of score and wts, where score is the propensity score values and 
#' wts are the corresponding stabilized inverse probability of treatment weights
#'
#' @export
#'
mvGPS <- function(D, C, common=FALSE, trim_w=FALSE, trim_quantile=0.99){
    m <- ncol(D)
    n <- nrow(D)
    if(common){
        if(is.list(C)) stop("common=TRUE, expecting C to be single matrix of common confounders", call.=FALSE)
        C <- rep(list(C), m)
    } else {
        if(!is.list(C)) stop("common=FALSE, C must be list of length m", call.=FALSE)
    }
    C_k <- unlist(lapply(C, ncol))
    C_n <- unlist(lapply(C, nrow))
    if(is.null(m)) stop("Exposure must be multivariate. See details to ensure formula is properly specified", call.=FALSE)
    if(!all(C_n==n)) stop("Each matrix in C must have same number of observations, n, as D", call.=FALSE)
    if(length(C)!=m) stop("Set of confounders not equal to number of exposures, m.")

    for(i in seq_len(m)){
        if(i==1){
            #marginal densities factorized
            d_1 <- lm(D[, i] ~ 1)
            d_1_mu <- coef(d_1)
            d_1_sigma <- summary(d_1)$sigma
            f_d_1 <- dnorm(D[, i], mean=d_1_mu, sd=d_1_sigma)

            gps_d_1 <- lm(D[, i] ~ C[[i]] - 1)
            gps_1_beta <- coef(gps_d_1)
            gps_1_Xb <- model.matrix(gps_d_1) %*% gps_1_beta
            gps_1_sigma <- summary(gps_d_1)$sigma
            f_gps_1 <- dnorm(D[, i], mean=gps_1_Xb, sd=gps_1_sigma)
        } else {
            cond_dens <- lapply(seq_len(m - 1) + 1, function(x){
                #full conditional marginal densities
                d_x <- lm(D[, x] ~ D[, seq_len(x-1)])
                d_x_beta <- coef(d_x)
                d_x_Xb <- model.matrix(d_x) %*% d_x_beta
                d_x_sigma <- summary(d_x)$sigma
                f_d_x <- dnorm(D[, x], mean=d_x_Xb, sd=d_x_sigma)

                #full conditional generalized propensity scores
                gps_x <- lm(D[, x] ~ D[, seq_len(x-1)] + C[[x]] - 1)
                gps_x_beta <- coef(gps_x)
                gps_x_Xb <- model.matrix(gps_x) %*% gps_x_beta
                gps_x_sigma <- summary(gps_x)$sigma
                f_gps_x <- dnorm(D[, x], gps_x_Xb, gps_x_sigma)

                return(list(marg=f_d_x, gps=f_gps_x))
            })
        }
    }
    cond_results <- unlist(cond_dens, recursive=FALSE)
    num_args <- cond_results[which(names(cond_results)=="marg")]
    num_args[["marg_1"]] <- f_d_1
    denom_args <- cond_results[which(names(cond_results)=="gps")]
    denom_args[["gps_1"]] <- f_gps_1
    
    score <- Reduce("*", denom_args)
    w <- Reduce("*", num_args)/score
    if(trim_w==TRUE){
        #trimming the large weights
        w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
        #trimming the small weights
        w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
    }
    return(list(score=score, wts=w))
}
