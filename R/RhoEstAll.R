#' @title RhoEstAll
#' @description
#'
#' @param yadd
#' @param ytime
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param h the bandwidth
#'
#' @return the correlation coefficient matrix

RhoEstAll <- function(yadd, ytime, u, u0, h, ktype = 'gaussian')  {
    p <- nrow(yadd)
    n1 <- length(u0)
    R <- matrix(0, nrow = p, ncol = n1)
    weight <- kernel_weight(u, u0 = u0, bw = h, poly_order = 0, ktype = ktype)
    for (i in 1:n1) {
        A <- yadd %*% weight[i, ]
        B <- ytime %*% weight[i, ]
        R[, i] <- Cubic(A, B)
    }
    return(R)
}
