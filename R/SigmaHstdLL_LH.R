#' @title SigmaHstdLL_LH
#'
#' @description  this routine computes the corresponding sigma matrix
#' with standard bandwidth at each u_i using local linear smoother.
#'
#' @param Y the observation matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @param h the bandwidth, scalar
#'
#' @return the estimator of  diagonal entries at each given u_i.
#' @export
#' @seealso \code{\link{CVLL}},\code{\link{CVHstdLL}}
#' @examples \dontrun{
#'
#' data(Yresid_ll)
#' data(u)
#' data(hstd_LL_LH)
#' # compute the estimator of diagonal entries
#' Sigma_LL_LH <- SigmaHstdLL_LH(Y = Yresid_ll, u = u, h = hstd_LL_LH$minimum)
#' Ystd_LL_LH = Yresid_ll / sqrt(Sigma_LL_LH)
#'}
SigmaHstdLL_LH  <- function(Y, u, h, ktype = 'gaussian')  {
    p <- nrow(Y)
    n <- ncol(Y)
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
    K <- kernelCompute(U_diff,  ktype = ktype, bw = h)
    UK_mat <- U_diff * K
    Y2 <- Y^2
    ExpU_mat <- exp(U_diff)
    CovDiag <- matrix(0, nrow = p, ncol = n)
    for (i in 1:n) {
        expu_i <- ExpU_mat[i, ]
        k_i <- K[i, ]
        uk_i <- UK_mat[i, ]
        SK <- sum(k_i)
        SUK <- sum(uk_i)
        v <- rep(0, p)
        for (j in 1:p) {
            y2_i <- Y2[j, ]
            v[j] <- Optimise_LMEBeta2(y2 = y2_i, expu = expu_i,
                                      k = k_i, uk = uk_i,
                                      SK = SK, SUK = SUK)
        }
        B <- Y2 / t(outer(expu_i, v, '^'))
        CovDiag[, i] <- B %*% k_i / SK
    }
    return(CovDiag)
}
