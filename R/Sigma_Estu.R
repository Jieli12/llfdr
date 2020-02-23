#' @title Sigma_Estu
#'
#' @param Y the observation matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @param h the bandwidth, scalar
#'
#' @return the diagonal entries estimator
#' @export
#'
#' @examples
Sigma_Estu <- function(Y, u, h, u0, ktype = 'gaussian') {
    weight <- kernel_weight(u, u0 = u0, bw = h, poly_order = 1, ktype = ktype)
    Y2 <- Y^2
    # sigma <- apply(weight, 1, .Product, Y = Y2)
    sigma <- apply(weight, 1, tcrossprod, y = Y2)
    sigma <- apply(sigma, 2, interpolate)
    return(sigma)
}
