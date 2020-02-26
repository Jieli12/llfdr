#' @title CvRhoAll
#' @description  This function computes the cross validation value of goal
#' function
#' @param yadd the sum
#' @param ytime the product
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param h the bandwidth
#' @export
#' @return the value of cross validation

CvRhoAll <- function(yadd, ytime, u, h, ktype = 'gaussian') {
    # p <- nrow(yadd)
    n <- ncol(yadd)
    # U <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    # U_diff <- U - t(U)
    U_diff <- computeUdiff(u)
    # U_diff_H <- U_diff / h
    # abs_diffh <- abs(U_diff_H)
    # K <- switch(ktype,
    #             gaussian = dnorm(U_diff_H) / h,
    #             epanech = {
    #                 ifelse(abs_diffh <= 1, 3 / 4 * ( 1 - U_diff_H^2) / h, 0)
    #             },
    #             triweight = {
    #                 ifelse(abs_diffh <= 1, 35 / 32 * ( 1 - U_diff_H^2)^3 / h, 0)
    #             },
    #             biweight = {
    #                 ifelse(abs_diffh <= 1, 15 / 16 * ( 1 - U_diff_H^2)^2 / h, 0)
    #             },
    #             tricube = {
    #                 ifelse(abs_diffh <= 1, 70 / 81 * ( 1 - U_diff_H^3)^3 / h, 0)
    #             },
    #             triangular = {
    #                 ifelse(abs_diffh <= 1, ( 1 - U_diff_H) / h, 0)
    #             },
    #             cosine = {
    #                 ifelse(abs_diffh <= 1, pi / 4 * cos( pi / 2 * U_diff_H) / h , 0)
    #             }
    # )
    K <- kernelCompute(U_diff, ktype = ktype, bw = h)
    SK <- rowSums(K) - diag(K)
    cv <- 0
    # R <- matrix(0, nrow = p, ncol = n)
    for(i in 1:n) {
        ai <- yadd[, i]
        bi <- ytime[, i]
        a_i <- yadd[, -i]
        b_i <- ytime[, -i]
        k_del_i <- K[i, -i]
        A <- a_i %*% k_del_i / SK[i]
        B <- b_i %*% k_del_i / SK[i]
        r <- Cubic(A, B)
        denom <- 1 - r^2
        cv <- cv + mean((ai - 2 * r * bi) / denom + log(denom))
    }
    return(cv)
}
