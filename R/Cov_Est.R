#' @title Cov_Est
#' @description this routine computes the estimator of covariance
#'
#' @param Y the p*n size matrix
#' @param weight the weight at u0
#'
#' @return the covariance matrix
#' @export
Cov_Est <- function(Y, weight) {
    # weight <- kernel_weight(u = u, u0 = u0, bw = h,
    #                         poly_order = poly_order, ktype = ktype)
    Cov <- xwxk1(Y, weight)
    # Cov <- tcrossprod(Y %*% diag(as.vector(weight)), Y)
    return(Cov)
}
