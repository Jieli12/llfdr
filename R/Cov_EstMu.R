#' @title Cov_EstMu
#' @description this routine computes the estimator of covariance
#'
#' @param Y the p*n size matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param h the bandwidth
#'
#' @return the covariance matrix
#' @export
Cov_EstMu <- function(Y, u, u0 = NULL, h, ktype = 'gaussian', poly_order = 0) {
    weight <- kernel_weight(u = u, u0 = u0, bw = h,
                            poly_order = poly_order, ktype = ktype)
    weight <- as.vector(weight)
    aver <- Y %*% weight
    Cov <- tcrossprod(Y %*% diag(weight), Y) - tcrossprod(aver)
    return(Cov)
}
