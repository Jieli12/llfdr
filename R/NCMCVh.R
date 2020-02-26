#' @title NCMCVh
#' @description This routine calculate the cross validation value using NCM method.
#' @param Y the observation, p * n matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param h the bandwidth, scalar
#' @export
#' @return the value of cross validation function

NCMCVh <- function(Y, u, h, ktype = 'gaussian') {
    p <- nrow(Y)
    n <- ncol(Y)
    weight <- kernel_weight(u, bw = h, poly_order = 0, ktype = ktype)
    cv <- 0
    for(i in 1:n) {
        covi <- CovDelete(Y[,-i], weight[i,-i])
        diff <- covi - tcrossprod(Y[,i])
        cv <- cv + norm(diff, "F")^2 / p
    }
    return(cv/n)
}
