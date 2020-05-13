#' Construct Covariate Balance Statistics for Simulation Models
#'
#' @param model_list character string identifying which methods to use when constructing weights. See details for a list of appropriate models
#' @param all_uni logical indicator. If TRUE then all univariate models specified
#' in model_list will be estimated for each exposure. If FALSE will only estimate weights
#' for the first exposure
#' @param D numeric matrix of dimension \code{n} by \code{2} designating values of the exposure
#' @param C numeric matrix of dimension \code{n} by \code{k} designating values of the confounders
#' @param trim_w logical indicator for whether to trim weights for mvGPS method. default is FALSE
#' 
#' @import WeightIt 
#' @import cobalt
#'
#' @export
#'
mvGPSsim_bal <- function(model_list, D, C, all_uni=TRUE, trim_w=FALSE){
    requireNamespace("WeightIt")
    requireNamespace("cobalt")
    k <- ncol(C)
    if(!ncol(D)>1) stop("'D' must be exposure matric with number of columns greater than or equal to 2", call.=FALSE)
    if(nrow(C)!=nrow(D)) stop("`C` and `D` must have same number of observations `n`")
    n <- nrow(D)
    model_list <- match.arg(model_list, c("mvGPS", "entropy", "CBPS", "PS", "GBM"), several.ok=TRUE)

    W <- lapply(model_list, function(mod){
       if(mod=="mvGPS"){
           w <- mvGPS(D~C, trim_w=trim_w)
           w_mvGPS <- list(w)
           names(w_mvGPS) <- "mvGPS"
           w <- w_mvGPS
       }else if(mod=="entropy"){
           if(all_uni==TRUE){
               w_entropy <- lapply(seq_len(ncol(D)), function(d){
                   w <- weightit(D[, d]~C, method="entropy")
                   w <- w$weights
               })
               names(w_entropy) <- paste0("entropy_D", seq_len(ncol(D)))
           } else {
               w_entropy <- weightit(D[, 1]~C, method="entropy")
               w_entropy <- list(w_entropy$weights)
               names(w_entropy) <- "entropy_D1"
           }
           w <- w_entropy
       }else if(mod=="CBPS"){
           if(all_uni==TRUE){
               w_CBPS <- lapply(seq_len(ncol(D)), function(d){
                   w <- weightit(D[, d]~C, method="cbps", over=FALSE)
                   w <- w$weights
               })
               names(w_CBPS) <- paste0("CBPS_D", seq_len(ncol(D)))
           } else {
               w_CBPS <- weightit(D[, 1]~C, method="cbps", over=FALSE)
               w_CBPS <- list(w_CBPS$weights)
               names(w_CBPS) <- "CBPS_D1"
           }
           w <- w_CBPS
       } else if(mod=="PS"){
           if(all_uni==TRUE){
               w_PS <- lapply(seq_len(ncol(D)), function(d){
                   w <- weightit(D[, d]~C, method="ps")
                   w <- w$weights
               })
               names(w_PS) <- paste0("PS_D", seq_len(ncol(D)))
           } else {
               w_PS <- weightit(D[, 1]~C, method="ps")
               w_PS <- list(w_PS$weights)
               names(w_PS) <- "PS_D1"
           }
           w <- w_PS
       } else if(mod=="GBM"){
           if(all_uni==TRUE){
               w_GBM <- lapply(seq_len(ncol(D)), function(d){
                   w <- weightit(D[, d]~C, method="GBM", stop.method="p.mean")
                   w <- w$weights
               })
               names(w_GBM) <- paste0("GBM_D", seq_len(ncol(D)))
           } else {
               w_GBM <- weightit(D[, 1]~C, method="GBM", stop.method="p.mean")
               w_GBM <- list(w_GBM$weights)
               names(w_GBM) <- "GBM_D1"
           }
           w <- w_GBM
       }
       return(w)
   })
    W <- unlist(W, recursive=FALSE, use.names=TRUE)

    #calculating correlation using weights
    D_corr_unweight <- lapply(seq_len(ncol(D)), function(dose){
        cobalt::col_w_corr(C, D[, dose])
    })
    names(D_corr_unweight) <- paste0("D", seq_len(ncol(D)), "_cor")
    D_corr_unweight_mat <- matrix(unlist(D_corr_unweight), nrow=k, ncol=ncol(D), dimnames=list(paste0("C", seq_len(k)), names(D_corr_unweight)))

    cor_list <- lapply(W, function(w){
        D_corr_w <- lapply(seq_len(ncol(D)), function(dose){
            cobalt::col_w_corr(C, D[, dose], weights=w)
        })
        names(D_corr_w) <- paste0("D", seq_len(ncol(D)), "_cor")
        D_corr_mat <- matrix(unlist(D_corr_w), nrow=k, ncol=ncol(D), dimnames=list(paste0("C", seq_len(k)), names(D_corr_w)))
        return(D_corr_mat)
    })
    cor_mat <- do.call(rbind, cor_list)
    cor_total <- rbind(D_corr_unweight_mat, cor_mat)
    cor_df <- data.frame(cor_total, method=rep(c("unweighted", names(cor_list)), each=k), covariate=paste0("C", seq_len(k)))

    #calculating summary balance metrics by group
    bal_group <- split(cor_df, cor_df$method)
    bal_metrics <- lapply(bal_group, function(x){
        euc_dist = sum(sqrt(sum(x[, grep("cor", names(x))]^2))) #since point of reference is [0, 0]
        max_cor = max(abs(x[, grep("cor", names(x))])) #this is equivalent to Chebyshev Distance
        avg_cor = mean(unlist(abs(x[, grep("cor", names(x))])))
        return(c(euc_dist, max_cor, avg_cor))
    })
    bal_metrics <- do.call(rbind, bal_metrics)
    colnames(bal_metrics) <- c("euc_dist", "max_cor", "avg_cor")
    bal_metrics <- data.frame(bal_metrics, method=names(bal_group))

    #calculating effective sample size by group
    ess <- lapply(W, function(w) sum(w)^2/sum(w^2))
    ess <- unlist(ess)

    return(list(W=W, cor_df=cor_df, bal_metrics=bal_metrics, ess=ess))
}
