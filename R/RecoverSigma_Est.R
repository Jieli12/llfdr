#' @title RecoverSigma_Est
#' @description This routine compute the i-th rows corresponding sigma for
#' revcoverring itself.
#'
#' @param Y the observation, p * n matrix
#' @param u the condtion, 1 * n matrix
#' @param h the bandwidth, scalar
#' @param u_test the center of u for test
#' @param ktype the kernel type, can be "gaussian", "epanech", "triweight",
#'              "biweight", "tricube", "triangular" and "cosine",
#'              the default of ktype is "gaussian".
#'
#' @return the  sigma

RecoverSigma_Est <- function(Y, u, h, u_test, ktype = 'gaussian') {
    u_test <- matrix(u_test, nrow = 1)
    w <- kernel_weight(u, u0 = u_test, bw = h, ktype = ktype)
    sigma <- t(tcrossprod(w, Y^2))
    return(sigma)
}
