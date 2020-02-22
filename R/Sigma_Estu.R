#' @title Sigma_Estu
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#' @param u_test the condition u for testing
#' @param ktype the kernel type, can be "gaussian", "epanech", "triweight",
#'              "biweight", "tricube", "triangular" and "cosine",
#'              the default of ktype is "gaussian".
#'
#' @return the diagonal entries estimator
#' @export
#'
#' @examples
Sigma_Estu <- function(Y, u, h, u_test, ktype = 'gaussian') {
    weight <- kernel_weight(u, u0 = u_test, bw = h, poly_order = 1, ktype = ktype)
    Y2 <- Y^2
    # sigma <- apply(weight, 1, .Product, Y = Y2)
    sigma <- apply(weight, 1, tcrossprod, y = Y2)
    sigma <- apply(sigma, 2, interpolate)
    return(sigma)
}
