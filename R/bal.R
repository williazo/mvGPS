#' Construct Covariate Balance Statistics for Models with Multivariate Exposure
#'
#'Assessing balance between exposure(s) and confounders is key when performing causal
#'analysis using propensity scores. We provide a list of several models to generate
#'weights to use in causal inference for multivariate exposures, and test the balancing property of these weights
#'using weighted Pearson correlations. In addition, returns the effective sample 
#'size.
#'
#' @inheritParams mvGPS
#' @param model_list character string identifying which methods to use when 
#' constructing weights. See details for a list of available models
#' @param all_uni logical indicator. If TRUE then all univariate models specified
#' in model_list will be estimated for each exposure. If FALSE will only estimate weights
#' for the first exposure
#' @param ... additional arguments to pass to \code{\link[WeightIt]{weightit}} function
#' if specifying one of these models in the model_list
#' 
#' @import WeightIt 
#' @import cobalt
#' @import gbm
#' @import CBPS
#' 
#' @details 
#' When using propensity score methods for causal inference it is crucial to 
#' check the balancing property of the covariates and exposure(s). To do this in
#' the multivariate case we first use a weight generating method from the available 
#' list shown below.
#' 
#' \subsection{Methods Available}{
#' 
#'   \itemize{
#'    \item "mvGPS": Multivariate generalized propensity score using Gaussian densities
#'    \item "entropy": Estimating weights using entropy loss 
#'    function without specifying propensity score \insertCite{tbbicke2020entropy}{mvGPS}
#'    \item "CBPS": Covariate balancing propensity score for continuous treatments
#'    which adds balance penalty while solving for propensity score parameters \insertCite{fong2018}{mvGPS}
#'    \item "PS": Generalized propensity score estimated using univariate Gaussian densities 
#'    \item "GBM": Gradient boosting to estimate the mean function of the propensity score, 
#'    but still maintains Gaussian distributional assumptions \insertCite{zhu_boosting}{mvGPS}
#'    }
#'    
#'    Note that only the \code{mvGPS} method is multivariate and all others are strictly univariate.
#'    For univariate methods weights are estimated for each exposure separately 
#'    using the \code{\link[WeightIt]{weightit}} function given the 
#'    confounders for that exposure in \code{C} when \code{all_uni=TRUE}. To estimate
#'    weights for only the first exposure set \code{all_uni=FALSE}.
#' }
#' 
#' It is also important to note that the weights for each method can be trimmed at 
#' the desired quantile by setting \code{trim_w=TRUE} and setting \code{trim_quantile}
#' in \[0.5, 1\]. Trimming is done at both the upper and lower bounds. For further details
#' see \code{\link[mvGPS]{mvGPS}} on how trimming is performed.
#' 
#' \subsection{Balance Metrics}{
#' 
#' In this package we include three key balancing metrics to summarize balance 
#' across all of the exposures.
#'  \itemize{
#'    \item Euclidean distance
#'    \item Maximum absolute correlation
#'    \item Average absolute correlation
#'   }
#'   \emph{Euclidean distance} is calculated using the origin point as reference, e.g. for \code{m=2} 
#'   exposures the reference point is \[0, 0\]. In this way we are calculating how far
#'   the observed set of correlation points are from perfect balance.
#'   
#'   \emph{Maximum absolute correlation} reports the largest single imbalance between
#'   the exposures and the set of confounders. It is often a key diagnostic as
#'   even a single confounder that is sufficiently out of balance can reduce performance.
#'   
#'   \emph{Average absolute correlation} is the sum of the exposure-confounder correlations.
#'   This metric summarizes how well, on average, the entire set of exposures is balanced.
#'  }
#'  
#'  \subsection{Effective Sample Size}{
#'  
#'  Effective sample size, ESS, is defined as
#'  \deqn{ESS=(\Sigma_i w_i)^{2}/\Sigma_i w_i^2,}
#'  where \eqn{w_i} are the estimated weights for a particular method \insertCite{kish_ess}{mvGPS}.
#'  Note that when \eqn{w=1} for all units that the \eqn{ESS} is equal to the sample size \eqn{n}.
#'  \eqn{ESS} decreases when there are extreme weights or high variability in the weights.
#'  }
#' 
#' @importFrom stats quantile
#' 
#' @examples 
#' #simulating data
#' sim_dt <- gen_D(method="u", n=150, rho_cond=0.2, s_d1_cond=2, s_d2_cond=2, 
#' k=3, C_mu=rep(0, 3), C_cov=0.1, C_var=1, d1_beta=c(0.5, 1, 0), 
#' d2_beta=c(0, 0.3, 0.75), seed=06112020)
#' D <- sim_dt$D
#' C <- sim_dt$C
#' 
#' #generating weights using mvGPS and potential univariate alternatives
#' require(WeightIt)
#' bal_sim <- bal(model_list=c("mvGPS", "entropy", "CBPS", "PS", "GBM"), D, 
#' C=list(C[, 1:2], C[, 2:3]))
#' 
#' #overall summary statistics
#' bal_sim$bal_metrics
#' 
#' #effective sample sizes
#' bal_sim$ess
#' 
#' #we can also trim weights for all methods 
#' bal_sim_trim <- bal(model_list=c("mvGPS", "entropy", "CBPS", "PS", "GBM"), D, 
#' C=list(C[, 1:2], C[, 2:3]), trim_w=TRUE, trim_quantile=0.9, p.mean=0.5)
#' #note that in this case we can also pass additional arguments using in 
#' #WeighIt package for entropy, CBPS, PS, and GBM such as specifying the p.mean
#' 
#' #can check to ensure all the weights have been properly trimmed at upper and 
#' #lower bound
#' all.equal(unname(unlist(lapply(bal_sim$W, quantile, 0.99))), 
#' unname(unlist(lapply(bal_sim_trim$W, max))))
#' all.equal(unname(unlist(lapply(bal_sim$W, quantile, 1-0.99))), 
#' unname(unlist(lapply(bal_sim_trim$W, min))))
#' 
#' @return 
#'    \itemize{
#'      \item \code{W}: list of weights generated for each model
#'      \item \code{cor_list}: list of weighted Pearson correlation coefficients for all confounders specified
#'      \item \code{bal_metrics}: data.frame with the Euclidean distance, maximum absolute correlation, and average absolute correlation by method
#'      \item \code{ess}: effective sample size for each of the methods used to generate weights
#'      \item \code{models}: vector of models used
#'      }
#'      
#' @references
#'     \insertAllCited{}
#'
#' @export
#'
bal <- function(model_list, D, C, common=FALSE, trim_w=FALSE, trim_quantile=0.99, all_uni=TRUE, ...){
    check_result <- D_C_check(D, C, common)
    assign("D", check_result$D)
    assign("C", check_result$C)
    
    m <- ncol(D)
    
    m_list <- as.list(
        match.arg(
            model_list,
            c("mvGPS", "entropy", "CBPS", "PS", "GBM"), 
            several.ok=TRUE
        )
    )
    if(any(is.na(pmatch(model_list, unlist(m_list))))){
        exclude_m <- model_list[is.na(pmatch(model_list, unlist(m_list)))]
        warning("The following model(s) specified but not used: ", paste(exclude_m, collapse=", "), call.=FALSE)
    }
    
    W <- lapply(m_list, function(mod){
       if(mod=="mvGPS"){
           w <- mvGPS(D, C, common = FALSE, trim_w, trim_quantile)
           #NOTE: we force common to be FALSE as we already create C from check in line 135
           w_mvGPS <- list(w$wts)
           names(w_mvGPS) <- "mvGPS"
           w <- w_mvGPS
       } else {
           #all other methods use the weightit function
           if(all_uni==TRUE){
               w <- lapply(seq_len(m), function(d){
                   w <- tryCatch(
                       expr= {
                           w <- weightit(D[, d]~C[[d]], method=mod, ...)
                           w <- w$weights
                           if(trim_w==TRUE){
                               #trimming the large weights
                               w <- ifelse(w<quantile(w, trim_quantile), w, quantile(w, trim_quantile))
                               #trimming the small weights
                               w <- ifelse(w>quantile(w, 1-trim_quantile), w, quantile(w, 1-trim_quantile))
                           }
                           return(w)
                       },
                       error = function(e) {
                           msg <- paste0(mod, " failed estimating exposure ", d, ". Returning null weights.")
                           warning(msg, call. = FALSE)
                           w <- NULL
                           return(w)
                       }
                   )
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
