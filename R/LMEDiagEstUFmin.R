#' @title LMEDiagEstUFmin
#'
#' @description  this routine computes the corresponding sigma matrix
#' with standard bandwidth at each u_i using local linear smoother.
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#' @param u0 the center condition, row matrix, size 1 * n2, the default is u.
#' @param ktype the kernel type, can be "gaussian", "epanech", "triweight",
#'              "biweight", "tricube", "triangular" and "cosine",
#'              the default of ktype is "gaussian".
#'
#' @return the estimator of  diagonal entries at each given u0.
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
LMEDiagEstUFmin  <- function(Y, u, u0, h, ktype = 'gaussian')  {
    p <- nrow(Y)
    n1 <- ncol(Y)
    n2 <- length(u0)
    u_rep <- rep(t(u), times = n2)
    u0_rep <- rep(t(u0), each = n1)
    U_diff <- matrix(u_rep - u0_rep, ncol = n1, byrow = TRUE)
    U_diff_H <- U_diff / h
    abs_diffh <- abs(U_diff_H)
    K <- switch(ktype,
                gaussian = dnorm(U_diff_H) / h,
                epanech = {
                    ifelse(abs_diffh <= 1, 3 / 4 * ( 1 - U_diff_H^2) / h, 0)
                },
                triweight = {
                    ifelse(abs_diffh <= 1, 35 / 32 * ( 1 - U_diff_H^2)^3 / h, 0)
                },
                biweight = {
                    ifelse(abs_diffh <= 1, 15 / 16 * ( 1 - U_diff_H^2)^2 / h, 0)
                },
                tricube = {
                    ifelse(abs_diffh <= 1, 70 / 81 * ( 1 - U_diff_H^3)^3 / h, 0)
                },
                triangular = {
                    ifelse(abs_diffh <= 1, ( 1 - U_diff_H) / h, 0)
                },
                cosine = {
                    ifelse(abs_diffh <= 1, pi / 4 * cos( pi / 2 * U_diff_H) / h , 0)
                }
    )
    UK_mat <- U_diff * K
    Y2 <- Y^2
    ExpU_mat <- exp(U_diff)
    CovDiag <- matrix(0, nrow = p, ncol = n)
    SK_rowsum <- rowSums(K)
    SUK_rowsum <- rowSums(UK_mat)
    for (i in 1:n2) {
        expu_i <- ExpU_mat[i, ]
        k_i <- K[i, ]
        uk_i <- UK_mat[i, ]
        SK <- SK_rowsum[i]
        SUK <- SUK_rowsum[i]
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
