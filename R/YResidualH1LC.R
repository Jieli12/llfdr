#' YResidualH1LC
#'
#' This function firstly estimates the local constant weighted
#' average at each center u(i) given bandwidth h. Here we use gaussian kerenel.
#' The return is the residual of Yi - \hat{Y}_{i}.
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#'
#' @return the residuals
#' @export
#'
#' @seealso \code{\link{YResidualH1LL}}
#' @examples \dontrun{
#'
#' data(Y)
#' data(u)
#' data(h1_lc)
#' Yresid_lc <- YResidualH1LC(Y = Y, u = u, h = h1_lc$minimum)
#' }
YResidualH1LC <- function(Y, u, h) {
    # p <- nrow(Y)
    # n <- ncol(Y)
    # yres <- matrix(0, nrow = p, ncol = n)
    # u_mat <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    # x <- (t(u_mat) - u_mat) / h
    # kernel <- dnorm(x) / h
    # sumk <- rowSums(kernel)
    # weight <- kernel / replicate(n, sumk)
    weight = kernel_weight(u, bw = h, poly_order = 0)
    # for (i in 1:n) {
    #     yres[,i] <- Y[,i] - Y %*% weight[i,]
    # }
    yres <- Y - apply(weight, 1, tcrossprod, y = Y)
    return(yres)
}
