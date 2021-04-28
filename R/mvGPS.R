#' Multivariate Generalized Propensity Score
#'
#' Estimate propensity scores for multivariate continuous exposure
#' by assuming joint normal conditional densities. Simultaneously estimate 
#' stabilized inverse probability of treatment weights (IPTW) using joint normal
#' density for marginal distribution of exposure.
#'
#' @param D numeric matrix of dimension \eqn{n} by \eqn{m} designating values of the exposures
#' @param C either a list of numeric matrices of length \eqn{m} of dimension 
#' \eqn{n} by \eqn{p_j} designating values of the confounders for each exposure 
#' value or if \code{common} is TRUE a single matrix of of dimension \eqn{n} by
#' \eqn{p} that represents common confounders for all exposures.
#' @param common logical indicator for whether C is a single matrix of common
#' confounders for all exposures. default is FALSE meaning C must be specified
#' as list of confounders of length \eqn{m}.
#' @param trim_w logical indicator for whether to trim weights. default is FALSE
#' @param trim_quantile numeric scalar used to specify the upper quantile to 
#' trim weights if applicable. default is 0.99
#' 
#' @details 
#' Generalized propensity scores (GPS) were proposed by 
#' \insertCite{hirano_continuous;textual}{mvGPS} and 
#' \insertCite{imai_causalGPS;textual}{mvGPS} to extend propensity scores to
#' handle continuous exposures. The GPS is constructed using the conditional 
#' density of the exposure given a set of confounders. These original methods 
#' and the subsequent literature have focused on the case of a single continuous
#' exposure where the GPS could be estimated using normal densities, kernel 
#' smoothing \insertCite{kennedy2017}{mvGPS}, gradient boosting 
#' \insertCite{zhu_boosting}{mvGPS}, and empirical likelihoods 
#' \insertCite{fong2018}{mvGPS}. In this package we provide an extension to this
#' literature to allow for multivariate continuous exposures to be estimated.
#' 
#' \subsection{Notation}{
#' 
#' Assume that we have a set of continuous exposures, \eqn{D}, of length 
#' \code{m}, i.e., \eqn{\mathbf{D}=D_{1}, \dots, D_{m}} collected on \eqn{n} 
#' units. Further, we assume that there exists a set of confounders 
#' \eqn{\mathbf{C}_{1},\dots,\mathbf{C}_{m}} for each
#' exposure of length \eqn{p_{j}} for \eqn{j=1,\dots,m}. The confounders are 
#' related to both the exposures and the outcome of interest such. Note that
#' the confounders need not be identical for all exposures. 
#' 
#' In our function we therefore expect that the argument \code{D} is a numeric
#' matrix of dimension \eqn{n\times m}, and that \code{C} is a list of length
#' \eqn{m} where each element is a matrix of dimension \eqn{n\times p_{j}}. For 
#' the case where we assume that all exposures have common confounders we 
#' set \code{common=TRUE} and \code{C} must be a matrix of dimension 
#' \eqn{n\times p}.
#' }
#' \subsection{mvGPS}{
#' 
#' We define the multivariate generalized propensity score, mvGPS, as
#' \deqn{mvGPS=f_{\mathbf{D}\mid \mathbf{C}_{1},\dots,\mathbf{C}_{m}}}
#' where \eqn{f} represents a joint multivariate conditional density function.
#' For our current development we specify \eqn{f} as multivariate normal, i.e.,
#' \deqn{\mathbf{D}\mid \mathbf{C}_{1},\dots,\mathbf{C}_{m}\sim N_{m}(\boldsymbol{\mu}, \Sigma).}
#' 
#' Factorizing this joint density we can rewrite the mvGPS as the product of 
#' univariate conditional densities, i.e.,
#' \deqn{mvGPS=f_{D_{m}\mid \mathbf{C}_{m}, D_{m-1},\dots,D_{1}}\times\cdots\times f_{D_{1}\mid\mathbf{C}_{1}}.}
#' We use this factorized version in our implementation, with parameters for each
#' distribution estimated through least squares.
#' }
#' \subsection{Constructing Weights}{
#' 
#' Following \insertCite{robins2000;textual}{mvGPS}, we use the mvGPS to 
#' construct stabilized inverse probability of treatment (IPTW) weights. These 
#' have been shown to balance confounders and return unbiased estimated of the
#' dose-response. Weights are constructed as
#' \deqn{w=\frac{f_{\mathbf{D}}}{f_{\mathbf{D}\mid \mathbf{C}_{1},\dots,\mathbf{C}_{m}}},}
#' where the marginal density \eqn{f_{\mathbf{D}}} of the exposures is assumed 
#' to be multivariate normal as well.
#' 
#' Writing the weights using completely factorized densities we have
#' \deqn{w=\frac{f_{D_{m}\mid D_{m-1},\dots, D_{1}}\times\cdots\times f_{D_{1}}}{f_{D_{m}\mid \mathbf{C}_{m}, D_{m-1},\dots,D_{1}}\times\cdots\times f_{D_{1}\mid\mathbf{C}_{1}}}.}
#' }
#' \subsection{Trimming}{
#' 
#' Often when using weights based on the propensity score, practitioners are 
#' concerned about the effect of extreme weights. It has been shown that an 
#' effective way to protect extreme weights is to trim them at a particular 
#' percentile \insertCite{crump2009dealing,lee2011weight}{mvGPS}. We allow
#' users to specify whether trimmed weights should be returned and at which
#' threshold. To trim weights set \code{trim_w=TRUE} and specify the desired
#' percentile as \code{trim_quantile=q}. Note that trimming is applied at 
#' both the upper and lower percentile thresholds, i.e.,
#' \deqn{w^{*}=w 1_{\{Q(w, 1-q)\le w \le Q(w, q)\}}+Q(w, 1-q) 1_{\{w < Q(w, 1-q)\}} + Q(w, q) 1_{\{w > Q(w, q)\}}}
#' }
#'
#' @importFrom stats model.frame model.matrix lm coef dnorm quantile
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' #generating confounded bivariate exposure
#' sim_dt <- gen_D(method="u", n=200, rho_cond=0.2, s_d1_cond=2, s_d2_cond=2, k=3, 
#' C_mu=rep(0, 3), C_cov=0.1, C_var=1, d1_beta=c(0.5, 1, 0), d2_beta=c(0, 0.3, 0.75), seed=06112020)
#' D <- sim_dt$D
#' C <- sim_dt$C
#' 
#' #generating weights and mvGPS
#' out_mvGPS <- mvGPS(D=D, C=list(C[, 1:2], C[, 2:3]))
#' 
#' # can apply trimming with default 99th percentile
#' out_mvGPS_trim <- mvGPS(D=D, C=list(C[, 1:2], C[, 2:3]), trim_w=TRUE)
#' 
#' #or assume all exposures have equivalent confounders
#' out_mvGPS_common <- mvGPS(D=D, C=C, common=TRUE)
#' 
#' @return list of score and wts, where score is the mvGPS score values and 
#' wts are the corresponding stabilized inverse probability of treatment weights
#'
#' @references
#'     \insertAllCited{}
#'
#' @export
#'
mvGPS <- function(D, C, common=FALSE, trim_w=FALSE, trim_quantile=0.99){
    check_result <- D_C_check(D, C, common)
    assign("D", check_result$D)
    assign("C", check_result$C)

    m <- ncol(D)
    
    for(i in seq_len(m)){
        if(i==1){
            #marginal densities factorized
            d_1 <- lm(D[, i] ~ 1)
            d_1_mu <- coef(d_1)
            d_1_sigma <- summary(d_1)$sigma
            f_d_1 <- dnorm(D[, i], mean=d_1_mu, sd=d_1_sigma)

            #generalized propensity score
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

#' Internal function for formatting and checking specification of exposures and
#' confounders
#' 
#' @inheritParams mvGPS
#' 
#' @keywords internal
#' 
#' @importFrom stats quantile
D_C_check <- function(D, C, common){
    D <- as.matrix(D)
    m <- ncol(D)
    if(is.null(colnames(D))) colnames(D) <- paste0("D", seq_len(m))
    n <- nrow(D)
    if(common){
        C <- as.matrix(C)
        if(is.null(colnames(C))) colnames(C) <- paste0("C", seq_len(ncol(C)))
        if(is.list(C)) stop("common=TRUE, expecting C to be single matrix of common confounders", call.=FALSE)
        C <- rep(list(C), m)
    } else {
        if(!is.list(C)) stop("common=FALSE, C must be list of length m", call.=FALSE)
        C <- lapply(seq_len(length(C)), function(x){
            X <- as.matrix(C[[x]])
             
            if(is.null(colnames(X))){
                p_j <- ncol(X)
                colnames(X) <- paste0("C_",x, "_", seq_len(p_j))
            }
            return(X)
            })
        
    }
    C_k <- unlist(lapply(C, ncol))
    C_n <- unlist(lapply(C, nrow))
    if(m < 2) stop("Exposure must be multivariate. See details to ensure formula is properly specified", call.=FALSE)
    if(!all(C_n==n)) stop("Each matrix in C must have same number of observations, n, as D", call.=FALSE)
    if(length(C)!=m) stop("Set of confounders not equal to number of exposures, m.")
    return(list(D=D, C=C))
}