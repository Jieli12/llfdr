#' Sigma_Gen_Model1
#' @description This runtine generates the p * p * length(u) size Sigma matrix.
#'
#' @param p the number of y variables
#' @param u the condition, it is a vector, one-dimensional array or
#' one-dimensional row(column) matrix
#'
#' @return p * p * length(u) size matrix
Sigma_Gen_Model1 <- function(p, u) {
    n <- length(u)
    Sigma <- array(0, dim = c(p, p, n))
    for (i in 1:p) {
        for (j in 1:p) {
            dnormu <- dnorm(u)
            Sigma[i, j, ] <- 0.5 * exp(u/2) *
                ((dnormu + 0.1) * as.numeric(abs(i - j)==1) + dnormu *
                     as.numeric(abs(i - j)==2) + 1 * as.numeric(i==j))
        }
    }
}
