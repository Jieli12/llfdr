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
#' @seealso \code{\link{CVLL}},\code{\link{CVHstdLL_LH}}
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

    U_diff <- computeUdiff(u)
    K <- kernelCompute(U_diff,  ktype = ktype, bw = h)
    UK_mat <- U_diff * K
    Y2 <- Y^2
    ExpU_mat <- exp(U_diff)
    CovDiag <- matrix(0, nrow = p, ncol = n)
    for (i in 1:n) {
        expu_i <- ExpU_mat[i, ]
        k_i <- K[i, ]
        uk_i <- UK_mat[i, ]
        mat_uk <- rbind(k_i, uk_i)
        SK <- sum(k_i)
        r0 <- sum(uk_i) / SK
        v <- rep(0, p)
        for (j in 1:p) {
            y2_i <- Y2[j, ]
            v[j] <- Optimise_LMEBeta2(y2 = y2_i, expu = expu_i,
                                      mat_uk = mat_uk, r0 = r0)
        }
        B <- Y2 / t(outer(expu_i, v, '^'))
        CovDiag[, i] <- B %*% k_i / SK
    }
    return(CovDiag)
}
