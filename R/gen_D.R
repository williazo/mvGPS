#' Generate Bivariate Multivariate Exposure
#'
#' Generate exposure from a bivariate normal distribution confounded by a set of 
#' variables \code{C}=\(C1, C2).
#'
#' @param method character value identifying which method to use when generating 
#' bivariate exposure. Options include "matrix_normal", "uni_cond", and "vector_normal".
#' See details for a brief explanation of each method. \code{uni_cond} is fastest
#' @param n integer value total number of units
#' @param rho_cond scalar value identifying conditional correlation of exposures given covariates between \[0, 1\]
#' @param s_d1_cond scalar value for conditional standard deviation of \code{D1}
#' @param s_d2_cond scalar value for conditional standard deviation of \code{D2}
#' @param k integer value determining number of covariates to generate in \code{C}.
#' @param C_mu numeric vector of mean values for covariates. Must be same length as \code{k}
#' @param C_cov scalar value representing constant correlation between covariates
#' @param C_var scalar value representing constant variance of covariates
#' @param d1_beta numeric vector of length \code{k} defining the mean of \code{D1} with respect to the covariates
#' @param d2_beta numeric vector of length \code{k} defining the mean of \code{D2} with respect to the covariates
#' @param seed integer value setting the seed of random generator to produce repeatable results. set to NULL by default
#'
#' @importFrom MASS mvrnorm
#' @importFrom matrixNormal rmatnorm J I vec
#' @importFrom stats rnorm cov2cor
#' 
#' @examples 
#' #generate bivariate exposures. D1 confounded by C1 and C2. D2 by C2 and C3
#' #uses univariate conditional normal to draw samples
#' sim_dt <- gen_D(method="u", n=200, rho_cond=0.2, s_d1_cond=2, s_d2_cond=2, k=3, 
#' C_mu=rep(0, 3), C_cov=0.1, C_var=1, d1_beta=c(0.5, 1, 0), d2_beta=c(0, 0.3, 0.75), seed=06112020)
#' D <- sim_dt$D
#' C <- sim_dt$C
#' 
#' #observed correlation should be close to true marginal value
#' cor(D); sim_dt$rho
#' 
#' 
#' #Use vector normal method instead of univariate method to draw samples
#' sim_dt <- gen_D(method="v", n=200, rho_cond=0.2, s_d1_cond=2, s_d2_cond=2, k=3, 
#' C_mu=rep(0, 3), C_cov=0.1, C_var=1, d1_beta=c(0.5, 1, 0), d2_beta=c(0, 0.3, 0.75), seed=06112020)
#' 
#' @details 
#'  \subsection{Generating Confounders}{
#'  
#'  We assume that there are a total of \code{k} confounders that are generated 
#'  from a multivariate normal distribution with equicorrelation covariance, i.e., 
#'  \deqn{\Sigma_{C}=\phi(\mathbf{1}\mathbf{1}^{T}-\mathbf{I})+\mathbf{I}\sigma^{2}_{C},}
#'  where \eqn{\mathbf{1}} is the column vector with all entries equal to 1,
#'  \eqn{\mathbf{I}} is the identity matrix, \eqn{\sigma^{2}_{C}} is a constant 
#'  standard deviation for all confounders, and \eqn{\phi} is the covariance of 
#'  any two confounders. Therefore, our random confounders 
#'  \code{C} follow the distribution
#'  \deqn{\mathbf{C}\sim N_{k}(\boldsymbol{\mu}_{C}, \Sigma_{C}).} 
#'  We draw a total of \code{n} samples from this multivariate normal distribution 
#'  using \code{\link[MASS]{mvrnorm}}.
#'  }
#'  \subsection{Generating Bivariate Exposure}{
#'  
#'  The first step when generating the bivariate exposure is to specify the 
#'  effects of the confounders \code{C}. We control this for each exposure value
#'  using the arguments \code{d1_beta} and \code{d2_beta} such that
#'  \deqn{E[D_{1}\mid \mathbf{C}]=\boldsymbol{\beta}^{T}_{D1}\mathbf{C}} and
#'  \deqn{E[D_{2}\mid \mathbf{C}]=\boldsymbol{\beta}^{T}_{D2}\mathbf{C}}.
#'  
#'  Note that by specifying \code{d1_beta} and \code{d2_beta} separately that the
#'  user can control the amount of overlap in the confounders for each exposure, 
#'  and how many of the variables in \code{C} are truly related to the exposures.
#'  For instance to have the exposure have identical confounding effects
#'  \code{d1_beta}=\code{d2_beta}, and they have separate confounding if there are 
#'  zero non-zero elements in common between \code{d1_beta} and \code{d2_beta}.
#'  
#'  To generate the bivariate conditional distribution of exposures given the set
#'  of confounders \code{C} we have the following three methods:
#'  \itemize{
#'   \item "matrix_normal"
#'   \item "uni_cond"
#'   \item "vector_normal"
#'   }
#'   
#'   "matrix_normal" uses the function \code{\link[matrixNormal]{rmatnorm}} to
#'   generate all \code{n} samples as
#'   \deqn{\mathbf{D}\mid\mathbf{C}\sim N_{n \times 2}(\boldsymbol{\beta}\mathbf{C}, \mathbf{I}_{n}, \Omega)}
#'   where \eqn{\boldsymbol{\beta}} is a column vector containing \eqn{\boldsymbol{\beta}^{T}_{D1}}
#'   and \eqn{\boldsymbol{\beta}^{T}_{D2}}, and \eqn{\Omega} is the conditional covariance matrix.
#'   
#'   "vector_normal" simply vectorizes the matrix_normal method above to generate
#'   a vector of length \eqn{n \times 2}. 
#'   
#'   "uni_cond" specifies the bivariate exposure using univariate conditional 
#'   factorization, which in the case of bivariate normal results in two univariate
#'   normal expressions.
#'   
#'   In general, we suggest using the univariate conditional, "uni_cond", method
#'   when generating exposures as it is substantially faster than both the 
#'   matrix normal and vector normal approaches.
#'   
#'   Note that the options use regular expression matching and can be specified 
#'   uniquely using either "m", "u", or "v".
#'  }
#'  \subsection{Marginal Covariance of Exposures}{
#'  
#'  As described above the exposures are drawn conditional on the set \code{C},
#'  so the marginal covariance of exposures is defined as 
#'  \deqn{\Sigma_{D}= \boldsymbol{\beta}\Sigma_{C}\boldsymbol{\beta}^{T}+\Omega.}
#'  In our function we return the true marginal covariance \eqn{\Sigma_{D}} as well
#'  as the true marginal correlation \eqn{\rho_{D}}.
#'  }
#'
#' @return 
#' \itemize{
#'  \item \code{D}: nx2 numeric matrix of the sample values for the exposures given the set \code{C}
#'  \item \code{C}: nxk numeric matrix of the sampled values for the confounding set \code{C}
#'  \item \code{D_Sigma}: 2x2 numeric matrix of the true marginal covariance of exposures
#'  \item \code{rho}: numeric scalar representing the true marginal correlation of exposures
#' }
#'
#' @export
gen_D <- function(method, n, rho_cond, s_d1_cond, s_d2_cond, k, C_mu, C_cov, C_var, d1_beta, d2_beta, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    method <- match.arg(method, c("matrix_normal", "uni_cond", "vector_normal"))

    C_Sigma <- I(k) * C_var + ((J(k) - I(k))*C_cov)
    C <- MASS::mvrnorm(n, mu=C_mu, Sigma=C_Sigma)
    colnames(C) <- paste0("C", seq_len(k))

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
        D <- matrix(D_vec, nrow=n, byrow=TRUE)
    }

    return(list(D=D, C=C, D_Sigma=D_Sigma, rho=rho))
}
