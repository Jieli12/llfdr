#' @title Cov_Est
#' @description this routine computes the estimator of covariance
#'
#' @param Y the p*n size matrix
#' @param u the condition
#' @param u_center the given center
#' @param h the bandwidth
#' @param ktype the kernel type, can be "gaussian", "epanech", "triweight",
#'              "biweight", "tricube", "triangular" and "cosine",
#'              the default of ktype is "gaussian".
#'
#' @return the covariance matrix
Cov_Est <- function(Y, u, u_center, h, ktype = 'gaussian') {
    weight <- kernel_weight(u = u, u0 = u_center, bw = h,
                            poly_order = 0, ktype = ktype)
    Cov <- tcrossprod(Y %*% diag(weight), Y)
    return(Cov)
}
