#' Construct Covariate Balance Statistics for Simulation Models
#'
#' @inheritParams mvGPS
#' @param model_list character string identifying which methods to use when constructing weights. See details for a list of appropriate models
#' @param all_uni logical indicator. If TRUE then all univariate models specified
#' in model_list will be estimated for each exposure. If FALSE will only estimate weights
#' for the first exposure
#' @param trim_w logical indicator for whether to trim weights for mvGPS method. default is FALSE
#' 
#' @import WeightIt 
#' @import cobalt
#' 
#' @importFrom stats quantile
#'
#' @export
#'
mvGPSsim_bal <- function(model_list, D, C, all_uni=TRUE, trim_w=FALSE, trim_quantile=0.99){
    requireNamespace("WeightIt")
    requireNamespace("cobalt")
    n <- nrow(D)
    m <- ncol(D)
    C_k <- unlist(lapply(C, ncol))
    C_n <- unlist(lapply(C, nrow))
    if(m<2) stop("'D' must be exposure matrix with number of columns, m, greater than or equal to 2", call.=FALSE)
    if(length(C)!=m) stop("`C` must be a list of length, m, where each list element represents the confounders for the corresponding exposure", call.=FALSE)
    if(!all(C_n==n)) stop("Each matrix in `C` must have same number of observations `n` as `D`", call.=FALSE)
    
    model_list <- match.arg(model_list, c("mvGPS", "entropy", "CBPS", "PS", "GBM"), several.ok=TRUE)

    W <- lapply(model_list, function(mod){
       if(mod=="mvGPS"){
           w <- mvGPS(D, C, trim_w=trim_w, trim_quantile=trim_quantile)
           w_mvGPS <- list(w)
           names(w_mvGPS) <- "mvGPS"
           w <- w_mvGPS
       }else if(mod=="entropy"){
           if(all_uni==TRUE){
               w_entropy <- lapply(seq_len(m), function(d){
                   w <- weightit(D[, d]~C[[d]], method="entropy")
                   w <- w$weights
                   if(trim_w==TRUE){
                       #trimming the large weights
                       w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                       #trimming the small weights
                       w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
                   }
                   w
               })
               names(w_entropy) <- paste0("entropy_D", seq_len(m))
           } else {
               w_entropy <- weightit(D[, 1]~C[[1]], method="entropy")
               w <- w_entropy$weights
               if(trim_w==TRUE){
                   #trimming the large weights
                   w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                   #trimming the small weights
                   w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
               }
               w_entropy <- list(w)
               names(w_entropy) <- "entropy_D1"
           }
           w <- w_entropy
       }else if(mod=="CBPS"){
           if(all_uni==TRUE){
               w_CBPS <- lapply(seq_len(m), function(d){
                   w <- weightit(D[, d]~C[[d]], method="cbps", over=FALSE)
                   w <- w$weights
                   if(trim_w==TRUE){
                       #trimming the large weights
                       w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                       #trimming the small weights
                       w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
                   }
                   w
               })
               names(w_CBPS) <- paste0("CBPS_D", seq_len(m))
           } else {
               w_CBPS <- weightit(D[, 1]~C[[1]], method="cbps", over=FALSE)
               w <- w_CBPS$weights
               if(trim_w==TRUE){
                   #trimming the large weights
                   w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                   #trimming the small weights
                   w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
               }
               w_CBPS <- list(w)
               names(w_CBPS) <- "CBPS_D1"
           }
           w <- w_CBPS
       } else if(mod=="PS"){
           if(all_uni==TRUE){
               w_PS <- lapply(seq_len(m), function(d){
                   w <- weightit(D[, d]~C[[d]], method="ps")
                   w <- w$weights
                   if(trim_w==TRUE){
                       #trimming the large weights
                       w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                       #trimming the small weights
                       w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
                   }
                   w
               })
               names(w_PS) <- paste0("PS_D", seq_len(m))
           } else {
               w_PS <- weightit(D[, 1]~C[[1]], method="ps")
               w <- w_PS$weights
               if(trim_w==TRUE){
                   #trimming the large weights
                   w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                   #trimming the small weights
                   w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
               }
               w_PS <- list(w)
               names(w_PS) <- "PS_D1"
           }
           w <- w_PS
       } else if(mod=="GBM"){
           if(all_uni==TRUE){
               w_GBM <- lapply(seq_len(m), function(d){
                   w <- weightit(D[, d]~C[[d]], method="GBM", stop.method="p.mean")
                   w <- w$weights
                   if(trim_w==TRUE){
                       #trimming the large weights
                       w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                       #trimming the small weights
                       w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
                   }
                   w
               })
               names(w_GBM) <- paste0("GBM_D", seq_len(m))
           } else {
               w_GBM <- weightit(D[, 1]~C[[1]], method="GBM", stop.method="p.mean")
               w <- w_GBM$weights
               if(trim_w==TRUE){
                   #trimming the large weights
                   w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                   #trimming the small weights
                   w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
               }
               w_GBM <- list(w)
               names(w_GBM) <- "GBM_D1"
           }
           w <- w_GBM
       }
       return(w)
   })
    W <- unlist(W, recursive=FALSE, use.names=TRUE)

    #calculating correlation using weights
    D_corr_unweight <- lapply(seq_len(m), function(dose){
        cobalt::col_w_corr(C[[dose]], D[, dose])
    })
    names(D_corr_unweight) <- paste0("D", seq_len(m), "_cor")
    # D_corr_unweight_mat <- matrix(unlist(D_corr_unweight), nrow=sum(C_k), ncol=m, 
    #                               dimnames=list(unlist(lapply(seq_len(m), function(i) paste0("C", i, "_", seq_len(C_k[i])))), names(D_corr_unweight)))

    cor_list <- lapply(W, function(w){
        D_corr_w <- lapply(seq_len(m), function(dose){
            cobalt::col_w_corr(C[[dose]], D[, dose], weights=w)
        })
        names(D_corr_w) <- paste0("D", seq_len(m), "_cor")
        # D_corr_mat <-  matrix(unlist(D_corr_w), nrow=sum(C_k), ncol=m, 
        #                       dimnames=list(c(paste0("C1_", seq_len(C_k[1])), paste0("C2_", seq_len(C_k[2]))), names(D_corr_w)))
        return(D_corr_w)
    })
    cor_list[["unweighted"]] <- D_corr_unweight
    # cor_df <- data.frame(cor_total, method=unlist(sapply(C_k, function(x) rep(c("unweighted", names(cor_list)), each=x))), 
    #                      covariate=c(paste0("C1_", seq_len(C_k[1])), paste0("C2_", seq_len(C_k[2]))))

    #calculating summary balance metrics by group
    # bal_group <- split(cor_df, cor_df$method)
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

    return(list(W=W, cor_list=cor_list, bal_metrics=bal_metrics, ess=ess))
}
