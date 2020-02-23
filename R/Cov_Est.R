#' @title Cov_Est
#' @description this routine computes the estimator of covariance
#'
#' @param Y the p*n size matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param h the bandwidth
#'
#' @return the covariance matrix
Cov_Est <- function(Y, u, u0 = NULL, h, ktype = 'gaussian', poly_order = 0) {
    weight <- kernel_weight(u = u, u0 = u0, bw = h,
                            poly_order = poly_order, ktype = ktype)
    Cov <- tcrossprod(Y %*% diag(weight), Y)
    return(Cov)
}
