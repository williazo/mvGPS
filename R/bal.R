#' Construct Covariate Balance Statistics for Models with Multivariate Exposure
#'
#' @inheritParams mvGPS
#' @param model_list character string identifying which methods to use when 
#' constructing weights. See details for a list of appropriate models
#' @param all_uni logical indicator. If TRUE then all univariate models specified
#' in model_list will be estimated for each exposure. If FALSE will only estimate weights
#' for the first exposure
#' @param trim_w logical indicator for whether to trim weights for mvGPS method. default is FALSE
#' 
#' @import WeightIt 
#' @import cobalt
#' 
#' @details 
#' When using propensity score methods for causal inference it is crucial to 
#' check the balancing property of the covariates and exposure(s). 
#' 
#' @importFrom stats quantile
#'
#' @export
#'
bal <- function(model_list, D, C, common=FALSE, trim_w=FALSE, trim_quantile=0.99, all_uni=TRUE){
    check_result <- D_C_check(D, C, common)
    e <- environment()
    list2env(check_result, e)
    
    m_list <- as.list(match.arg(model_list, 
                                    c("mvGPS", "entropy", "CBPS", "PS", "GBM"), 
                                    several.ok=TRUE))
    if(any(is.na(pmatch(model_list, unlist(m_list))))){
        exclude_m <- model_list[is.na(pmatch(model_list, unlist(m_list)))]
        warning("The following model(s) specified but not used: ", paste(exclude_m, collapse=", "), call.=FALSE)
    }
    
    W <- lapply(m_list, function(mod){
       if(mod=="mvGPS"){
           w <- mvGPS(D, C, common, trim_w, trim_quantile)
           w_mvGPS <- list(w$wts)
           names(w_mvGPS) <- "mvGPS"
           w <- w_mvGPS
       } else { #all other methods use the weightit function
           if(all_uni==TRUE){
               w <- lapply(seq_len(m), function(d){
                   w <- weightit(D[, d]~C[[d]], method=mod)
                   w <- w$weights
                   if(trim_w==TRUE){
                       #trimming the large weights
                       w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                       #trimming the small weights
                       w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
                   }
                   w
               })
               names(w) <- paste0(mod, "_", colnames(D))
           } else {
               w <- weightit(D[, 1]~C[[1]], method=mod)
               w <- w$weights
               if(trim_w==TRUE){
                   #trimming the large weights
                   w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                   #trimming the small weights
                   w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
               }
               w <- list(w)
               names(w) <- paste0(mod, "_", colnames(D)[1])
           }
       }
       return(w)
   })
    W <- unlist(W, recursive=FALSE, use.names=TRUE)

    #calculating correlation using weights
    D_corr_unweight <- lapply(seq_len(m), function(dose){
        col_w_corr(C[[dose]], D[, dose])
    })
    names(D_corr_unweight) <- paste0(colnames(D), "_cor")
    
    cor_list <- lapply(W, function(w){
        D_corr_w <- lapply(seq_len(m), function(dose){
            col_w_corr(C[[dose]], D[, dose], weights=w)
        })
        names(D_corr_w) <- paste0(colnames(D), "_cor")
        return(D_corr_w)
    })
    cor_list[["unweighted"]] <- D_corr_unweight

    #calculating summary balance metrics by group
    bal_metrics <- lapply(cor_list, function(x){
        cor_vec <- unlist(x)
        euc_dist = sum(sqrt(sum(cor_vec^2))) #since point of reference is [0, 0]
        max_cor = max(abs(cor_vec)) #this is equivalent to Chebyshev Distance
        avg_cor = mean(unlist(abs(cor_vec)))
        return(c(euc_dist, max_cor, avg_cor))
    })
    bal_metrics <- do.call(rbind, bal_metrics)
    colnames(bal_metrics) <- c("euc_dist", "max_cor", "avg_cor")
    bal_metrics <- data.frame(bal_metrics, method=row.names(bal_metrics))

    #calculating effective sample size by group
    ess <- lapply(W, function(w) sum(w)^2/sum(w^2))
    ess <- unlist(ess)

    return(list(W=W, cor_list=cor_list, bal_metrics=bal_metrics, ess=ess, 
                models=unlist(m_list)))
}
