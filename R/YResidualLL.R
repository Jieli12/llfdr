#' YResidualLL
#'
#' This function firstly estimates the local linear weighted
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
YResidualLL  <- function(Y, u, h) {
    p <- nrow(Y)
    n <- ncol(Y)
    yres <- matrix(0, nrow = p, ncol = n)

    u_mat <- replicate(n, u)
    u_diff <- t(u_mat) - u_mat
    x <- u_diff / h
    kernel <- dnorm(x) / h
    u_diff_2 <- u_diff^2
    sn1 <- rowSums(kernel * u_diff)
    sn2 <- rowSums(kernel * u_diff_2)
    Sn1 <- replicate(n, sn1)
    Sn2 <- replicate(n, sn2)
    Weight <- kernel * (Sn2 - u_diff * Sn1)
    Weight_sum = rowSums(Weight)
    Weight = Weight / replicate(n, Weight_sum)
    for (i in 1:n) {
        yres[,i] <- Y[,i] - Y %*% Weight[i,]
    }
    return(yres)
}
