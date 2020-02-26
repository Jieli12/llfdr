#' @param YY the data, num_non_zero * n matrix
#' @param Cov_vec the weighted average of YY
#' @param u the condition
#' @param u_test the test center
#' @param h the bandwidth
#' @param ysquare the square of standadization of y
#' @param Sigmau the diagonal entries estimator at u_test
#' @param Cov_LL the covariance matrix
#' @param ktype the kernel type, can be "gaussian", "epanech", "triweight",
#'              "biweight", "tricube", "triangular" and "cosine",
#'              the default of ktype is "gaussian".
#'
#' @title ShrinkageCoefLC
#' @description  This function computes the shrinkage covariance after threshold .
#' @export
#' @return the shrinkage covariance


ShrinkageCoefLC <- function(YY, Cov_vec, u, u_test, h, ysquare, Sigmau, Cov_LL,
                      ktype = 'gaussian') {
    w <- kernel_weight(u, u0 = u_test, bw = h, poly_order = 0, ktype = ktype)
    n <- length(w)
    #  Shrinkage
    #  Compute the value of \hat{m}(u)
    meanvar <- mean(diag(Sigmau))
    prior <- diag(meanvar / diag(Sigmau))

    # Compute the numerator of coefficient, i.e., \hat{\beta}^2(u), here we denote
    # it as phi
    phi <- 0
    for(i in 1:n) {
        diff_offdiagonal <- w[i] * (YY[, i] - Cov_vec)
        diff_diagonal <- w[i] * (ysquare[, i] - 1)
        phi <- phi + sum(diff_offdiagonal^2) * 2 + sum(diff_diagonal^2)
    }
    # Compute the denominator of coefficient, i.e., \hat{\alpha}^2(u), here we
    # denote it as gamma.
    gamma <- 2 * sum(Cov_vec^2) + sum((1 - diag(prior))^2)

    # Compute the shrinkage covariance.
    coef <- phi / (gamma + phi)
    Cov1 <- coef * prior + (1 - coef) * Cov_LL
    return(Cov1)
}

