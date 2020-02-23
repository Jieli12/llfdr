#' @title RecoverSigma_Est
#' @description This routine compute the i-th rows corresponding sigma for
#' revcoverring itself.
#'
#' @param Y the observation, p * n matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param h the bandwidth, scalar
#' @param ktype the kernel type, can be "gaussian", "epanech", "triweight",
#'              "biweight", "tricube", "triangular" and "cosine",
#'              the default of ktype is "gaussian".
#'
#' @return the  sigma

RecoverSigma_Est <- function(Y, u, h, u0, ktype = 'gaussian') {
    # u0 <- matrix(u0, nrow = 1)
    w <- kernel_weight(u, u0 = u0, bw = h, ktype = ktype)
    sigma <- t(tcrossprod(w, Y^2))
    return(sigma)
}
