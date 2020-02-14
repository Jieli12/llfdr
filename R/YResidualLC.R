#' YResidualLC
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
#' @examples
YResidualLC <- function(Y, u, h) {
    p <- nrow(Y)
    n <- ncol(Y)
    yres <- matrix(0, nrow = p, ncol = n)
    u_mat <- replicate(n, u)
    x <- (t(u_mat) - u_mat) / h
    kernel <- dnorm(x) / h
    sumk <- rowSums(kernel)
    weight <- kernel / replicate(n, sumk)
    for (i in 1:n) {
        yres[,i] <- Y[,i] - Y %*% weight[i,]
    }
    return(yres)
}
