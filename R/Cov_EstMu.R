#' @title Cov_EstMu
#' @description this routine computes the estimator of covariance
#'
#' @param Y the p*n size matrix
#' @param weight the bandwidth at u0
#'
#' @return the covariance matrix
#' @export
Cov_EstMu <- function(Y, weight) {
    weight <- as.vector(weight)
    aver <- Y %*% weight
    # Cov <- tcrossprod(Y %*% diag(weight), Y) - tcrossprod(aver)
    Cov <- xwxk1(Y, weight) - tcrossprod(aver)
    return(Cov)
}
