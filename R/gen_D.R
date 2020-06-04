#' Generate Random Multivariate Exposure
#'
#' This function allows the user to generate exposure from a multivariate normal distribution confounded by a set of variables \code{C}.
#'
#' @param method character value identifying which method so use when generating exposure. for bivariate model \code{uni_cond} is fastest
#' @param n integer value total number of units
#' @param rho_cond scalar value identifying conditional correlation of exposures given covariates between \[0, 1\]
#' @param s_d1_cond scalar value for conditional standard deviation of \code{D1}
#' @param s_d2_cond scalar value for conditional standard deviation of \code{D2}
#' @param num_covs integer value determining number of covariates to generate
#' @param C_mu numeric vector of mean values for covariates. must be same length as \code{num_covs}
#' @param C_cov scalar value representing constant correlation between covariates
#' @param C_var scalar value indicating constant variance of covariates
#' @param d1_beta numeric vector of length \code{num_covs} defining the mean of \code{D1} with respect to the covariates
#' @param d2_beta numeric vector of length \code{num_covs} defining the mean of \code{D2} with respect to the covariates
#' @param seed integer value setting the seed of random generator to produce repeatable results. set to NULL by default
#'
#' @importFrom MASS mvrnorm
#' @importFrom matrixNormal rmatnorm J I vec
#' @importFrom stats rnorm cov2cor
#'
#' @export
gen_D <- function(method, n, rho_cond, s_d1_cond, s_d2_cond, num_covs, C_mu, C_cov, C_var, d1_beta, d2_beta, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    method <- match.arg(method, c("matrix_normal", "uni_cond", "vector_normal"))
    requireNamespace("MASS")
    requireNamespace("matrixNormal")

    C_Sigma <- I(num_covs) * C_var + ((J(num_covs) - I(num_covs))*C_cov)
    C <- MASS::mvrnorm(n, mu=C_mu, Sigma=C_Sigma)
    colnames(C) <- paste0("C", seq_len(num_covs))

    d_beta <- cbind(d1_beta, d2_beta)
    d_xbeta <- C %*% d_beta
    d1_xbeta <- d_xbeta[, 1]
    d2_xbeta <- d_xbeta[, 2]

    s_d_cond <- c(s_d1_cond, s_d2_cond)

    D_corr_cond <- I(2) + matrix(c(0, rho_cond, rho_cond, 0), nrow=2, ncol=2, byrow=TRUE)
    D_Sigma_cond <- outer(s_d_cond, s_d_cond) * D_corr_cond
    
    # calculating the true marginal covariance
    D_Sigma <- t(d_beta)%*%C_Sigma%*%d_beta + D_Sigma_cond
    # true marginal correlation matrix
    D_corr <- cov2cor(D_Sigma)
    # true marginal correlation
    rho <- D_corr[1, 2]

    if(method=="matrix_normal"){
        D <- rmatnorm(M=d_xbeta, U=I(n), V=D_Sigma_cond)
    } else if(method=="uni_cond"){
        D1 <- rnorm(n, d1_xbeta, s_d1_cond)
        D2 <- rnorm(n, d2_xbeta + (s_d2_cond/s_d1_cond)*rho_cond*(D1-d1_xbeta), sqrt((1-rho_cond^2)*s_d2_cond^2))
        D <- cbind(D1, D2)
    } else if(method=="vector_normal"){
        D_vec <- MASS::mvrnorm(1, mu=vec(d_xbeta), Sigma=kronecker(I(n), D_Sigma_cond))
        D <- matrix(D_vec, nrow=n)
    }

    return(list(D=D, C=C, D_Sigma=D_Sigma, rho=rho))
}
