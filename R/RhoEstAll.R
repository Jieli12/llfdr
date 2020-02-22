#' @title RhoEstAll
#' @description
#'
#' @param yadd
#' @param ytime
#' @param u
#' @param u_test
#' @param h
#' @param ktype
#'
#' @return the correlation coefficient matrix

RhoEstAll <- function(yadd, ytime, u, u_test, h, ktype = 'gaussian')  {
    p <- nrow(yadd)
    n1 <- length(u_test)
    R <- matrix(0, nrow = p, ncol = n1)
    weight = kernel_weight(u, u0 = u_test, bw = h, poly_order = 0, ktype = ktype)
    for (i in 1:n1) {
        A <- yadd %*% weight[i, ]
        B <- ytime %*% weight[i, ]
        R[, i] <- Cubic(A, B)
    }
    return(R)
}
