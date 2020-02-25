#' YResidualH1LL
#'
#' This function firstly estimates the local linear weighted
#' average at each center u(i) given bandwidth h. Here we use gaussian kerenel.
#' The return is the residual of Yi - \\hat{Y}_{i}.
#'
#' @param Y the observation matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @param h the bandwidth, scalar
#'
#' @return the residuals
#' @export
#'
#' @seealso \code{\link{YResidualH1LC}}
#' @examples \dontrun{
#'
#' data(Y)
#' data(u)
#' data(h1_ll)
#' Yresid_ll <- YResidualH1LL(Y = Y, u = u, h = h1_ll$minimum)
#' }
YResidualH1LL  <- function(Y, u, h, ktype = 'gaussian') {
    # p <- nrow(Y)
    # n <- ncol(Y)
    # yres <- matrix(0, nrow = p, ncol = n)

    # u_mat <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    # u_diff <- t(u_mat) - u_mat
    # x <- u_diff / h
    # kernel <- dnorm(x) / h
    # u_diff_2 <- u_diff^2
    # sn1 <- rowSums(kernel * u_diff)
    # sn2 <- rowSums(kernel * u_diff_2)
    # Sn1 <- replicate(n, sn1)
    # Sn2 <- replicate(n, sn2)
    # weight <- kernel * (Sn2 - u_diff * Sn1)
    # weight_sum = rowSums(weight)
    # weight = weight / replicate(n, weight_sum)
    weight <- kernel_weight(u, bw = h, ktype = ktype, poly_order = 1)
    # for (i in 1:n) {
    #     yres[,i] <- Y[,i] - Y %*% weight[i,]
    # }
    yres <- Y - apply(weight, 1, tcrossprod, y = Y)
    return(yres)
}
