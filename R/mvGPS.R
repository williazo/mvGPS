#' Multivariate Generalized Propensity Score (mvGPS)
#'
#' Construct weights using propensity scores for multivariate continuous exposure
#' by assuming joint normal marginal and conditional densities. See details for
#' description of weights and additional methodology details.
#'
#' @param formula an object of class formula that describes the generalized
#' propensity score model. The outcome must be contain at least two exposures
#' @param data data.frame that contains variables used in formula
#' @param trim_w logical indicator for whether to trim weights. default is FALSE
#' @param trim_quantile numeric scalar used to specify the upper quantile to 
#' trim weights if applicabile. default is 0.99
#'
#' @importFrom stats model.frame model.matrix lm coef dnorm
#'
#' @export
#'
mvGPS <- function(formula, data = sys.frame(sys.parent()), trim_w=FALSE, trim_quantile=0.99){
    if(class(formula) != "formula"){
        stop("`formula` must be a formula", call. = FALSE)
    } else if (length(formula) != 3){
        stop("`formula` must be 2-sided", call. = FALSE)
    }

    y <- model.frame(formula, data)[, 1]
    X <- model.matrix(formula, data)
    if(ncol(y) < 2) stop("exposure must be multivariate. see details to ensure formula is properly specified")

    for(i in seq_len(ncol(y))){
        if(i==1){
            #marginal densities factorized
            d_1 <- lm(y[, i] ~ 1)
            d_1_mu <- coef(d_1)
            d_1_sigma <- summary(d_1)$sigma
            f_d_1 <- dnorm(y[, i], mean=d_1_mu, sd=d_1_sigma)

            gps_d_1 <- lm(y[, i] ~ X - 1)
            gps_1_beta <- coef(gps_d_1)
            gps_1_Xb <- X %*% gps_1_beta
            gps_1_sigma <- summary(gps_d_1)$sigma
            f_gps_1 <- dnorm(y[, i], mean=gps_1_Xb, sd=gps_1_sigma)
        } else {
            cond_dens <- lapply(seq_len(ncol(y)-1)+1, function(x){
                #full conditional marginal densities
                d_x <- lm(y[, x] ~ y[, seq_len(x-1)])
                d_x_beta <- coef(d_x)
                d_x_Xb <- model.matrix(d_x) %*% d_x_beta
                d_x_sigma <- summary(d_x)$sigma
                f_d_x <- dnorm(y[, x], mean=d_x_Xb, sd=d_x_sigma)

                #full conditional generalized propensity scores
                gps_x <- lm(y[, x] ~ y[, seq_len(x-1)] + X - 1)
                gps_x_beta <- coef(gps_x)
                gps_x_Xb <- model.matrix(gps_x) %*% gps_x_beta
                gps_x_sigma <- summary(gps_x)$sigma
                f_gps_x <- dnorm(y[, x], gps_x_Xb, gps_x_sigma)

                return(list(marg=f_d_x, gps=f_gps_x))
            })
        }
    }
    cond_results <- unlist(cond_dens, recursive=FALSE)
    num_args <- cond_results[which(names(cond_results)=="marg")]
    num_args[["marg_1"]] <- f_d_1
    denom_args <- cond_results[which(names(cond_results)=="gps")]
    denom_args[["gps_1"]] <- f_gps_1

    Reduce("*", num_args)

    w <- Reduce("*", num_args)/Reduce("*", denom_args)
    if(trim_w==TRUE){
        #trimming the large weights
        w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
        #trimming the small weights
        w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
    }
    return(w)
}
