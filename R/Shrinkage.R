#' @title Shrinkage
#' @description  This function computes the shrinkage covariance after threshold .
#'
#' @param Y the data, p * n matrix
#' @param w the weight
#' @param Cov the covariance when the i-th deleted
#' @param Ind_t0 the indicator matrix
#'
#' @return the shrinkage covariance
Shrinkage <- function(Y, w, Cov, Ind_t0) {
    p <- nrow(Y)
    n <- ncol(Y)

    #  Shrinkage
    #  Compute the value of \hat{m}(u)
    meanvar <- mean(diag(Cov))
    prior <- meanvar * diag(p)

    # Compute the numerator of coefficient, i.e., \hat{\beta}^2(u), here we denote
    # it as phi
    phi <- 0
    for(i in 1:n) {
        diff <- w[i] * ( (tcrossprod(Y[, i]) - Cov) * Ind_t0)
        phi <- phi + norm(diff, "F")^2
    }
    # Compute the denominator of coefficient, i.e., \hat{\alpha}^2(u), here we
    # denote it as gamma.
    gamma <- norm(Cov - prior,"F")^2

    # Compute the shrinkage covariance.
    coef <- phi / (gamma + phi)
    Cov1 <- coef * prior + (1 - coef) * Cov
    return(Cov1)
}
